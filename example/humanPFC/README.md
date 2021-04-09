# Human prefrontal cortex snm3C-seq analysis
This is an example of analyzing 4238 cells from human prefrontal cortex. It includes [embedding](## Clustering) and compartment calling at 100kb resolution, domain calling at 25kb resolution, loop calling at 10kb resolution.
## Prepare directory
```bash
for r in 10k 25k 100k; do for c in `seq 1 22`; do mkdir -p cell_matrix/${r}b_resolution/chr${c}; mkdir -p imputed_matrix/${r}b_resolution/chr${c}; done; mkdir imputed_matrix/${r}b_resolution/merged/; done
```
## Generate matrix at multiple resolutions


## Clustering
### Impute at 100kb resolution
```bash
# parallelize at cell level
cell=$(cut -f1 cell_4238_meta_cluster.txt | sed '1d' | head -${SGE_TASK_ID} | tail -1)
res0=100k
res=100000
for c in `seq 1 22`; 
do 
command time hicluster impute-cell --indir cell_matrix/${res0}b_resolution/chr${c}/ --outdir imputed_matrix/${res0}b_resolution/chr${c}/ --cell ${cell} --chrom ${c} --res ${res} --pad 1 --chrom_file hg19.autosomal.chrom.sizes; 
done
```
### Generate embedding
```bash
for c in `seq 1 22`; do awk -v c=$c '{printf("imputed_matrix/100kb_resolution/chr%s/%s_chr%s_pad1_std1_rp0.5_sqrtvc.hdf5\n",c,$1,c)}' cell_4238_meta_cluster.txt > imputed_matrix/100kb_resolution/filelist/imputelist_pad1_std1_rp0.5_sqrtvc_chr${c}.txt; echo $c; done
# parallelize at chromosome level
c=${SGE_TASK_ID}
command time hicluster embed-concatcell-chr --cell_list imputed_matrix/100kb_resolution/filelist/imputelist_pad1_std1_rp0.5_sqrtvc_chr${c}.txt --outprefix imputed_matrix/100kb_resolution/merged/embed/pad1_std1_rp0.5_sqrtvc_chr${c} --res ${res}

# merge chromosomes together
ls imputed_matrix/100kb_resolution/merged/embed/*npy > imputed_matrix/100kb_resolution/filelist/embedlist_pad1_std1_rp0.5_sqrtvc.txt
command time hicluster embed-mergechr --embed_list imputed_matrix/100kb_resolution/filelist/embedlist_pad1_std1_rp0.5_sqrtvc.txt --outprefix imputed_matrix/100kb_resolution/merged/embed/pad1_std1_rp0.5_sqrtvc
```
### Plot result

## Loop calling
### Impute at 10kb resolution
We start by impute matrices at 10kb resolution. The running time of random walk increases cubically with the dimension of a matrix, which makes the imputation at 10kb resolution very slow for large chromosomes. However, at 10kb resolution, we usually only care about the contacts within certain distance (e.g. < 10 Mb). This allowed us to use a sliding window method to speed up the computation. Specifically, we only compute random walk within the square matrix of dimension w, and move downstreamly by step size of s. The imputation results of these windows are merged to generate the final imputation. Empirically when computing a matrix > 12000 dimensions (human chr1-chr12), the sliding window of size 30-40M will significantly accelerate the imputation and returns accurate results of contacts within 10M.  
After imputation, the matrix is normalized by the decay with a Z-score transform at each insertion distance, and the differences between each pixel and its local background is also saved for later t-tests of loop candidates.
```bash
# parallelize at cell level
cell=$(cut -f1 cell_4238_meta_cluster.txt | sed '1d' | head -${SGE_TASK_ID} | tail -1)
res0=10k
res=10000
for c in `seq 1 12`; 
do 
command time hicluster impute-cell --indir cell_matrix/${res0}b_resolution/chr${c}/ --outdir imputed_matrix/${res0}b_resolution/chr${c}/ --cell ${cell} --chrom ${c} --res ${res} --chrom_file hg19.autosomal.chrom.sizes --pad 2 --output_dist 10050000 --window_size 30000000 --step_size 10000000 --mode pad2_std1_rp0.5_sqrtvc; 
command time hicluster loop-bkg-cell --indir imputed_matrix/${res0}b_resolution/ --cell ${cell} --chrom ${c} --res ${res} --impute_mode pad2_std1_rp0.5_sqrtvc; 
done
for c in `seq 13 22`; 
do 
command time hicluster impute-cell --indir cell_matrix/${res0}b_resolution/chr${c}/ --outdir imputed_matrix/${res0}b_resolution/chr${c}/ --cell ${cell} --chrom ${c} --res ${res} --chrom_file hg19.autosomal.chrom.sizes --pad 2 --output_dist 10050000 --mode pad2_std1_rp0.5_sqrtvc; 
command time hicluster loop-bkg-cell --indir imputed_matrix/${res0}b_resolution/ --cell ${cell} --chrom ${c} --res ${res} --impute_mode pad2_std1_rp0.5_sqrtvc; 
done
```
### Sum cell
The next step is to sum the matrices of single cells in the same group.
```

```

### write cell list
```bash
for i in `seq 0 4`; do for c in `seq 1 7`; do for a in p q; do awk -v c=$c$a '{printf("/gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/chr%s/%s_chr%s_pad2_std1_rp0.5_sqrtvc\n",c,$1,c)}' celllist_covgroup${i}.txt > L23_covgroup${i}_pad2_std1_rp0.5_sqrtvc_chr${c}${a}_looplist.txt; done; done; done
for i in `seq 0 4`; do for c in `seq 1 22`; do awk -v c=$c '{printf("/gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/chr%s/%s_chr%s_pad2_std1_rp0.5_sqrtvc\n",c,$1,c)}' celllist_covgroup${i}.txt > L23_covgroup${i}_pad2_std1_rp0.5_sqrtvc_chr${c}_looplist.txt; done; done
```
### merge cells by coverage groups
```bash
cat <(for w in ws30 ws40; do for c in `seq 1 12`; do for i in `seq 0 4`; do echo $i $c pad2_std1_rp0.5_${w}; done; done; done) | sort -k2,2n -k1,1rn > paralist_ws.txt 
cat <(for c in `seq 1 22`; do for i in `seq 0 4`; do echo $i $c pad2_std1_rp0.5_sqrtvc; done; done) | sort -k2,2n -k1,1rn > paralist.txt 
cat <(for c in `seq 1 7`; do for a in p q; do for i in `seq 0 4`; do echo $i ${c}${a} pad2_std1_rp0.5_sqrtvc; done; done; done) | sort -k2,2 -k1,1rn > paralist_split.txt 
```
### write group list
```bash
for c in `seq 1 22`; do ls /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup?_pad2_std1_rp0.5_sqrtvc_dist_trim/L23_covgroup?_pad2_std1_rp0.5_sqrtvc_dist_trim_chr${c}.hdf5 | sed 's/.hdf5//g' > filelist/L23_pad2_std1_rp0.5_sqrtvc_chr${c}_grouplist.txt; done
for c in `seq 1 12`; do for w in 30 40; do ls /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup?_pad2_std1_rp0.5_ws${w}_dist_trim/L23_covgroup?_pad2_std1_rp0.5_ws${w}_dist_trim_chr${c}.hdf5 | sed 's/.hdf5//g' > filelist/L23_pad2_std1_rp0.5_ws${w}_chr${c}_grouplist.txt; done; done
for c in `seq 1 7`; do for a in p q; do ls /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup?_pad2_std1_rp0.5_sqrtvc_dist_trim/L23_covgroup?_pad2_std1_rp0.5_sqrtvc_dist_trim_chr${c}${a}.hdf5 | sed 's/.hdf5//g' > filelist/L23_pad2_std1_rp0.5_sqrtvc_chr${c}${a}_grouplist.txt; done; done
```
### merge groups
```bash
cat <(cat paralist.txt paralist_ws.txt | cut -d' ' -f2,3 | sort -k1,1n -k2,2 -u) <(cat paralist_split.txt | cut -d' ' -f2,3 | sort -k1,1 -k2,2 -u) > paralist_group.txt

file=paralist_group.txt
c=$(cut -f1 -d' ' ${file} | head -${SGE_TASK_ID} | tail -1)
mode=$(cut -f2 -d' ' ${file} | head -${SGE_TASK_ID} | tail -1)
command time python /gale/ddn/snm3C/humanPFC/code/loop_sumcell_chr.py --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/filelist/L23_${mode}_chr${c}_looplist.txt --group_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/filelist/L23_${mode}_chr${c}_grouplist.txt --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_${mode}_dist_trim/L23_${mode}_dist_trim_chr${c} --res 10000
```

