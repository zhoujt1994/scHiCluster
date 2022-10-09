import pathlib

import numpy as np
import pandas as pd
import ray
import xarray as xr


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


@ray.remote
def save_sample_chunk(
        output_path,
        cool_ds_paths,
        loop_ds_path,
        sample_chunk,
        chrom,
        da_name,
        value_types,
):
    ds = load_cool_ds_chrom(cool_ds_paths, chrom)
    loop_ds = xr.open_zarr(f"{loop_ds_path}/{chrom}/")
    loop_mask = loop_ds["loop"] > 0

    bin1_chunk_size, bin2_chunk_size, *_ = ds[da_name].encoding["chunks"]
    bin1_idx = ds.get_index("bin1")
    bin2_idx = ds.get_index("bin2")

    first = True
    cur_loops = 0
    for bin1_start in range(0, bin1_idx.size, bin1_chunk_size):
        bin1_chunk = bin1_idx[bin1_start: bin1_start + bin1_chunk_size]
        ds_chunks = []
        for bin2_start in range(0, bin2_idx.size, bin2_chunk_size):
            bin2_chunk = bin2_idx[bin2_start: bin2_start + bin2_chunk_size]
            loop_chunk = loop_mask.sel(
                {"bin1": bin1_chunk, "bin2": bin2_chunk}
            ).to_pandas()
            loop_x, loop_y = np.where(loop_chunk)
            loop_chrom_global_x = loop_x + bin1_start
            loop_chrom_global_y = loop_y + bin2_start
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

            # save zarr
            loop_index = pd.RangeIndex(cur_loops, cur_loops + loop_x.size, name="loop")
            da = xr.DataArray(loop_pixel_data,
                              dims=["loop", "sample_id", f"{da_name}_value_type"],
                              coords={'sample_id': sample_chunk,
                                      f"{da_name}_value_type": value_types,
                                      "loop": loop_index})
            da.encoding["chunks"] = (50000, sample_chunk.size, 1)
            chunk_ds = xr.Dataset({da_name: da,
                                   'loop_bin1_id': pd.Series(loop_chrom_global_x, index=loop_index),
                                   'loop_bin2_id': pd.Series(loop_chrom_global_y, index=loop_index)})
            cur_loops += loop_x.size
            ds_chunks.append(chunk_ds)
        if len(ds_chunks) > 0:
            ds_to_save = xr.concat(ds_chunks, dim="loop")
            if first:
                ds_to_save.to_zarr(output_path, mode="w")
                first = False
            else:
                ds_to_save.to_zarr(output_path, append_dim="loop")
    return


def extract_loop_ds(
        cool_ds_paths, loop_ds_path, output_path, da_name, chroms, value_types
):
    pathlib.Path(output_path).mkdir(exist_ok=True, parents=True)
    value_types = pd.Index(value_types)

    futures = []
    for chrom in chroms:
        ds = load_cool_ds_chrom(cool_ds_paths, chrom)
        sample_idx = ds.get_index("sample_id")
        *_, sample_chunk_size, _ = ds[da_name].encoding["chunks"]

        for sample_start in range(0, sample_idx.size, sample_chunk_size):
            sample_chunk = sample_idx[sample_start: sample_start + sample_chunk_size]
            chunk_path = f'{output_path}/{chrom}_chunk_{sample_start}'
            future = save_sample_chunk.remote(
                output_path=chunk_path,
                cool_ds_paths=cool_ds_paths,
                loop_ds_path=loop_ds_path,
                sample_chunk=sample_chunk,
                chrom=chrom,
                da_name=da_name,
                value_types=value_types
            )
            futures.append(future)
    ray.get(futures)
    return
