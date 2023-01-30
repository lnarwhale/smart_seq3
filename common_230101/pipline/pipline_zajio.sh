source bulk.sh

work_dir=$(pwd)
cd ..
all_dir=$(pwd)
cd ${work_dir}
pair="single"
species="mouse"
group_1_name="X1C"
group_2_name="X2C"
group1="X1C"
group2="X2C"
sample="${group1} ${group2}"
reference_rna_dir="/data1/shenluemou/reference/mapping/hisat2/mouse/grcm38/mouse_grcm38_rRNA"
reference_dir="/data1/shenluemou/reference/mapping/hisat2/mouse/grcm38/mouse_grcm38"
gtf_dir="/data1/shenluemou/reference/gtf/mouse/grcm38/mouse_grcm38.gtf"
bed_dir="/data1/shenluemou/reference/bed/mouse/grcm38/mouse_grcm38.bed"
#----------cut--------------------
mkdir ../result/
mkdir ../result/log
cd ../result/log
log=$(pwd)
cd $work_dir
mkdir ../result/cut
for name in ${sample};do
	smartseq3_cut ${all_dir}/data ${name} ${all_dir}/result/cut
	echo "the cut of "${name}" finished"
done
#-------fastp---------------
mkdir ../result/fastp
for name in ${sample};do
	fastp_control ${all_dir}/result/cut ${name} ${all_dir}/result/fastp
	echo "the fastp of ${name} finished"
done
#---------trim----------------
mkdir ../result/trim
for name in ${sample};do
	trim_g ${all_dir}/result/fastp ${name} ${all_dir}/result/trim 50
	echo "the trim_galore of "${name}" finished"
done

#---------qurrna---------------
mkdir ../result/1de
for name in ${sample};do
	de_rna ${all_dir}/result/trim  ${name} ${all_dir}/result/1de ${reference_rna_dir} 
	echo "the de_rrna of "${name}" finished"
done

#-------bw--------------------
mkdir ../result/2map
for name in ${sample};do
	fastqtobw ${all_dir}/result/1de ${name} ${reference_dir} ${all_dir}/result/2map
	echo "the mapping of "${name}" finished"
done

#------cager_all--------------
mkdir ../result/2_1cager
cp ${all_dir}/result/2map/*_sorted.bam ${all_dir}/result/2_1cager
cager_a ${all_dir}/result/2_1cager ${species} ${all_dir}/result/2_1cager

#--------------rnaqc-------------
mkdir ../result/3qc
for name in ${sample};do
	rnaqc ${all_dir}/result/2map ${name} ${all_dir}/result/3qc ${gtf_dir} ${bed_dir}
	echo "the bam-qc of "${name}" finished"
done

#---------------abundance----------
mkdir ../result/4abundance
for name in ${sample};do
	estimate_abundance ${all_dir}/result/2map ${name} ${gtf_dir} ${all_dir}/result/4abundance
	echo "the abundance of "${name}" finished"
done
mkdir ../result/5matrix
total_abundance ${all_dir}/result/4abundance ${all_dir}/result/5matrix
read_matrix ${all_dir}/result/5matrix ${all_dir}/result/4abundance

#----------------saturation------------
saturation ${all_dir}/result/5matrix/qualimeta.txt ${all_dir}/result/5matrix/

#----------------PCA------------------
touch ${all_dir}/result/5matrix/matrix.csv
echo "Sample,Group" > ${all_dir}/result/5matrix/matrix.csv
for name in ${group1};do
        echo "${group_1_name},${name}" >> ${all_dir}/result/5matrix/matrix.csv
done
for name in ${group2};do
        echo "${group_2_name},${name}" >> ${all_dir}/result/5matrix/matrix.csv
done
Rscript pca_plot.R -m ${all_dir}/result/5matrix/matrix.csv -g ${all_dir}/result/5matrix/gene_fpkm_matrix.csv -o ${all_dir}/result/5matrix


