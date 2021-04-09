# Human prefrontal cortex snm3C-seq analysis

# write cell list

for i in `seq 0 4`; do for c in `seq 1 7`; do for a in p q; do awk -v c=$c$a '{printf("/gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/chr%s/%s_chr%s_pad2_std1_rp0.5_sqrtvc\n",c,$1,c)}' celllist_covgroup${i}.txt > L23_covgroup${i}_pad2_std1_rp0.5_sqrtvc_chr${c}${a}_looplist.txt; done; done; done
for i in `seq 0 4`; do for c in `seq 1 22`; do awk -v c=$c '{printf("/gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/chr%s/%s_chr%s_pad2_std1_rp0.5_sqrtvc\n",c,$1,c)}' celllist_covgroup${i}.txt > L23_covgroup${i}_pad2_std1_rp0.5_sqrtvc_chr${c}_looplist.txt; done; done

# merge cells by coverage groups

cat <(for w in ws30 ws40; do for c in `seq 1 12`; do for i in `seq 0 4`; do echo $i $c pad2_std1_rp0.5_${w}; done; done; done) | sort -k2,2n -k1,1rn > paralist_ws.txt 
cat <(for c in `seq 1 22`; do for i in `seq 0 4`; do echo $i $c pad2_std1_rp0.5_sqrtvc; done; done) | sort -k2,2n -k1,1rn > paralist.txt 
cat <(for c in `seq 1 7`; do for a in p q; do for i in `seq 0 4`; do echo $i ${c}${a} pad2_std1_rp0.5_sqrtvc; done; done; done) | sort -k2,2 -k1,1rn > paralist_split.txt 

# write group list

for c in `seq 1 22`; do ls /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup?_pad2_std1_rp0.5_sqrtvc_dist_trim/L23_covgroup?_pad2_std1_rp0.5_sqrtvc_dist_trim_chr${c}.hdf5 | sed 's/.hdf5//g' > filelist/L23_pad2_std1_rp0.5_sqrtvc_chr${c}_grouplist.txt; done
for c in `seq 1 12`; do for w in 30 40; do ls /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup?_pad2_std1_rp0.5_ws${w}_dist_trim/L23_covgroup?_pad2_std1_rp0.5_ws${w}_dist_trim_chr${c}.hdf5 | sed 's/.hdf5//g' > filelist/L23_pad2_std1_rp0.5_ws${w}_chr${c}_grouplist.txt; done; done
for c in `seq 1 7`; do for a in p q; do ls /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup?_pad2_std1_rp0.5_sqrtvc_dist_trim/L23_covgroup?_pad2_std1_rp0.5_sqrtvc_dist_trim_chr${c}${a}.hdf5 | sed 's/.hdf5//g' > filelist/L23_pad2_std1_rp0.5_sqrtvc_chr${c}${a}_grouplist.txt; done; done

# merge groups

cat <(cat paralist.txt paralist_ws.txt | cut -d' ' -f2,3 | sort -k1,1n -k2,2 -u) <(cat paralist_split.txt | cut -d' ' -f2,3 | sort -k1,1 -k2,2 -u) > paralist_group.txt

file=paralist_group.txt
c=$(cut -f1 -d' ' ${file} | head -${SGE_TASK_ID} | tail -1)
mode=$(cut -f2 -d' ' ${file} | head -${SGE_TASK_ID} | tail -1)
command time python /gale/ddn/snm3C/humanPFC/code/loop_sumcell_chr.py --cell_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/filelist/L23_${mode}_chr${c}_looplist.txt --group_list /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/filelist/L23_${mode}_chr${c}_grouplist.txt --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_${mode}_dist_trim/L23_${mode}_dist_trim_chr${c} --res 10000


mode=pad2_std1_rp0.5_sqrtvc
cd merged/L23_${mode}_dist_trim/
for c in `seq 13 22`; do ln -s /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_pad2_std1_rp0.5_sqrtvc_dist_trim/L23_pad2_std1_rp0.5_sqrtvc_dist_trim_chr${c}.loop.hdf5 L23_${mode}_dist_trim_chr${c}.loop.hdf5; done

command time python /gale/ddn/snm3C/humanPFC/code/loop_mergechr.py --inprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_${mode}_dist_trim/L23_${mode}_dist_trim --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_${mode}_dist_trim/L23_${mode}_dist_trim --chrom_file /gale/netapp/home/zhoujt/genome/hg19/hg19.autosomal.chrom.sizes

command time python /gale/ddn/snm3C/humanPFC/code/loop_mergechr.py --inprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_${mode}_dist_trim/L23_${mode}_dist_trim --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_${mode}_dist_trim/L23_${mode}_dist_trim.split --chrom_file /gale/netapp/home/zhoujt/genome/hg19/hg19.autosomal.chrom.sizes --split_file /gale/netapp/home/zhoujt/genome/hg19/hg19.chrsplit.bed

for i in range(5):
	for c in chrom:
		loop = pd.read_hdf(f'/gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup{i}_pad2_std1_rp0.5_sqrtvc_dist_trim/L23_covgroup{i}_pad2_std1_rp0.5_sqrtvc_dist_trim_chr{c}.loop.hdf5', key='chr'+c)
		loop.to_hdf(f'/gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup{i}_pad2_std1_rp0.5_sqrtvc_dist_trim/L23_covgroup{i}_pad2_std1_rp0.5_sqrtvc_dist_trim_chr{c}.loop.hdf5', key='loop', mode='w')
		print(i, c)

mode=pad2_std1_rp0.5_sqrtvc
for i in `seq 0 4`; do command time python /gale/ddn/snm3C/humanPFC/code/loop_mergechr.py --inprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup${i}_${mode}_dist_trim/L23_covgroup${i}_${mode}_dist_trim --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup${i}_${mode}_dist_trim/L23_covgroup${i}_${mode}_dist_trim --chrom_file /gale/netapp/home/zhoujt/genome/hg19/hg19.autosomal.chrom.sizes; done
for i in `seq 0 4`; do command time python /gale/ddn/snm3C/humanPFC/code/loop_mergechr.py --inprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup${i}_${mode}_dist_trim/L23_covgroup${i}_${mode}_dist_trim --outprefix /gale/ddn/snm3C/humanPFC/smoothed_matrix/10kb_resolution/merged/L23_covgroup${i}_${mode}_dist_trim/L23_covgroup${i}_${mode}_dist_trim.split --chrom_file /gale/netapp/home/zhoujt/genome/hg19/hg19.autosomal.chrom.sizes --split_file /gale/netapp/home/zhoujt/genome/hg19/hg19.chrsplit.bed; done

