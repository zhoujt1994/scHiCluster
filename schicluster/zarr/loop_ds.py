import pathlib

import dask.array as da
import numpy as np
import pandas as pd
import ray
import xarray as xr
import zarr


def load_cool_ds_chrom(paths, chrom):
    chrom_ds_paths = [f"{path}/{chrom}" for path in paths]
    ds = xr.open_mfdataset(
        chrom_ds_paths,
        engine="zarr",
        concat_dim="sample_id",
        combine="nested",
        decode_cf=False,
    )
    return ds


def init_empty_loop_array(
        cool_ds, loop_mask, value_types, da_name, output_path, chrom, loop_chunk_size=50000
):
    bin1_chunk_size, bin2_chunk_size, sample_chunk_size, _ = cool_ds[da_name].encoding[
        "chunks"
    ]
    bin1_idx = cool_ds.get_index("bin1")
    bin2_idx = cool_ds.get_index("bin2")
    sample_idx = cool_ds.get_index("sample_id")

    loop_x = []
    loop_y = []
    # get loop coords, and make sure the order is by bin chunk,
    # the same order as when iterate bin and save real data
    for bin1_start in range(0, bin1_idx.size, bin1_chunk_size):
        bin1_chunk = bin1_idx[bin1_start: bin1_start + bin1_chunk_size]
        for bin2_start in range(0, bin2_idx.size, bin2_chunk_size):
            bin2_chunk = bin2_idx[bin2_start: bin2_start + bin2_chunk_size]
            loop_chunk = loop_mask.sel(
                {"bin1": bin1_chunk, "bin2": bin2_chunk}
            ).to_pandas()
            _loop_x, _loop_y = np.where(loop_chunk)
            if _loop_x.size > 0:
                loop_x.append(_loop_x + bin1_start)
                loop_y.append(_loop_y + bin2_start)

    loop_x = np.concatenate(loop_x)
    loop_y = np.concatenate(loop_y)
    n_loop = loop_x.size

    empty_loop_array = xr.DataArray(
        da.zeros(
            (n_loop, sample_idx.size, value_types.size),
            chunks=(loop_chunk_size, sample_chunk_size, 1),
            dtype="float32",
        ),
        dims=["loop", "sample_id", f"{da_name}_value_type"],
        coords={"sample_id": sample_idx, f"{da_name}_value_type": value_types},
    )
    loop_idx = empty_loop_array.get_index("loop")
    empty_loop_ds = xr.Dataset({da_name: empty_loop_array})
    empty_loop_ds.coords['loop_bin1_id'] = pd.Series(loop_x, index=loop_idx)
    empty_loop_ds.coords['loop_bin2_id'] = pd.Series(loop_y, index=loop_idx)
    empty_loop_ds.to_zarr(f"{output_path}/{chrom}", mode="w")
    return


@ray.remote
def save_sample_chunk(
        output_path,
        cool_ds_paths,
        loop_position_ds_path,
        sample_chunk,
        sample_start,
        chrom,
        da_name,
        value_types,
):
    ds = load_cool_ds_chrom(cool_ds_paths, chrom)
    loop_ds = xr.open_zarr(f"{loop_position_ds_path}/{chrom}/")
    loop_mask = loop_ds["loop"] > 0
    output_da = zarr.open(f"{output_path}/{chrom}/{da_name}")

    bin1_chunk_size, bin2_chunk_size, *_ = ds[da_name].encoding["chunks"]
    bin1_idx = ds.get_index("bin1")
    bin2_idx = ds.get_index("bin2")

    print(f'Saving {chrom} Sample {sample_start}-{sample_start + sample_chunk.size}')

    cur_loop_idx_start = 0
    for bin1_start in range(0, bin1_idx.size, bin1_chunk_size):
        bin1_chunk = bin1_idx[bin1_start: bin1_start + bin1_chunk_size]
        for bin2_start in range(0, bin2_idx.size, bin2_chunk_size):
            bin2_chunk = bin2_idx[bin2_start: bin2_start + bin2_chunk_size]
            loop_chunk = loop_mask.sel(
                {"bin1": bin1_chunk, "bin2": bin2_chunk}
            ).to_pandas()
            loop_x, loop_y = np.where(loop_chunk)
            if loop_x.size == 0:
                continue

            # turn 2D matrix shape into 1D pixel shape
            # matrix_data.shape = [bin1, bin2, sample, value_type]
            matrix_data = (
                ds[da_name]
                .sel(
                    {
                        "bin1": bin1_chunk,
                        "bin2": bin2_chunk,
                        "sample_id": sample_chunk,
                        f"{da_name}_value_type": value_types,
                    }
                )
                .values
            )
            # loop_pixel_data.shape = [loop, sample, value_type]
            loop_pixel_data = matrix_data[loop_x, loop_y, ...]

            # save data to zarr at chunk location
            cur_loop_idx_end = cur_loop_idx_start + loop_x.size
            sample_end = sample_start + sample_chunk.size
            output_da[
            cur_loop_idx_start:cur_loop_idx_end, sample_start:sample_end, :
            ] = loop_pixel_data

            # update loop idx
            cur_loop_idx_start = cur_loop_idx_end
    return


def create_loop_ds(
        cool_ds_paths,
        loop_position_ds_path,
        output_path,
        da_name,
        chroms,
        value_types,
        min_loop_count=1,
        loop_chunk_size=50000
):
    pathlib.Path(output_path).mkdir(exist_ok=True, parents=True)
    value_types = pd.Index(value_types)

    futures = []
    for chrom in chroms:
        print(f'Init {chrom} empty zarr dataset')
        ds = load_cool_ds_chrom(cool_ds_paths, chrom)
        loop_ds = xr.open_zarr(f"{loop_position_ds_path}/{chrom}/")
        loop_mask = loop_ds["loop"] >= min_loop_count
        init_empty_loop_array(
            cool_ds=ds,
            loop_mask=loop_mask,
            value_types=value_types,
            da_name=da_name,
            output_path=output_path,
            chrom=chrom,
            loop_chunk_size=loop_chunk_size
        )

        # save sample chunks in parallel
        sample_idx = ds.get_index("sample_id")
        *_, sample_chunk_size, _ = ds[da_name].encoding["chunks"]
        for sample_start in range(0, sample_idx.size, sample_chunk_size):
            sample_chunk = sample_idx[sample_start: sample_start + sample_chunk_size]
            # noinspection PyArgumentList
            future = save_sample_chunk.remote(
                output_path=output_path,
                cool_ds_paths=cool_ds_paths,
                loop_position_ds_path=loop_position_ds_path,
                sample_chunk=sample_chunk,
                sample_start=sample_start,
                chrom=chrom,
                da_name=da_name,
                value_types=value_types,
            )
            futures.append(future)
    ray.get(futures)
    return
