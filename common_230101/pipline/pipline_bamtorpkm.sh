source bulk.sh

work_dir=$(pwd)
cd ..
all_dir=$(pwd)
species="mouse"
sample="nss_za_X0C nss_za_X2C smart_mou_X0C smart_mou_X2C_3 smart_za_X0C smart_za_X2C"
gtf_dir="/data1/shenluemou/reference/gtf/mouse/grcm38/mouse_grcm38.gtf"
bed_dir="/data1/shenluemou/reference/bed/mouse/grcm38/mouse_grcm38.bed"

cd $work_dir
#--------------rnaqc-------------------
mkdir ../result_bamtorpkm
mkdir ../result_bamtorpkm/1qc
for name in ${sample};do
	rnaqc ${all_dir}/bam $name ${all_dir}/result_bamtorpkm/1qc $gtf_dir $bed_dir
	echo "the bam-qc of "$name" has been finished"
done

#--------------abundance---------------
mkdir ../result_bamtorpkm/2abundance
for name in ${sample};do
	estimate_abundance ${all_dir}/bam $name $gtf_dir $all_dir/result_bamtorpkm/2abundance
	echo "the abundance of "${name}" finished"
done
mkdir ../result_bamtorpkm/3matrix
total_abundance $all_dir/result_bamtorpkm/2abundance $all_dir/result_bamtorpkm/3matrix
read_matrix $all_dir/result_bamtorpkm/3matrix $all_dir/result_bamtorpkm/2abundance







 
