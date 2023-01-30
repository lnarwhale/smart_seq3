#!/bin/bash
function help_bulk_or_pair(){
	echo "bop_fastp_control | 目的:对数据进行fastp质控｜位置:1.fastq所在文件夹 2.fastq的名字 3.质控后数据的
输出文件夹"
	echo "------------------------------------"
	echo "bop_map | 目的:对数据进行比对|位置:1.fastq_dir 2.name 3.reference_dir&name(index) 4.out_dir"
	echo "------------------------------------"
	echo "" 
}
#---------------------------------------------
function bop_fastp_control(){
	#--function:fastp quality control
	#--input:1.fastq_dir 2.name 3.out_dir
	#--need:fastp
	fastq_dir=$1
	name=$2
	out_dir=$3
	echo "the fastp of "$fastq_dir/$name" begin"
	/data1/shenluemou/biosoft/fastp \
	-i $fastq_dir/$name"_R1.fq" \
	-I $fastq_dir/$name"_R2.fq" \
	-o $out_dir/$name"_R1.fq" \
	-O $out_dir/$name"_R2.fq" \
	--html $out_dir/$name".html" \
	--trim_poly_g --poly_g_min_len 5 \
	--trim_poly_x --poly_x_min_len 5 \
	--cut_front --cut_tail --cut_window_size 4 \
	--qualified_quality_phred 15 \
	--low_complexity_filter \
	--complexity_threshold 30 \
	--length_required 30 \
	--thread 4
}
#----------------------------------------------
function bop_map(){
	#--function:from fastq to map to get sam and tran-to bam and bw
	#--input:1.fastq_dir 2.name 3.reference_dir&name(index) 4.out_dir
	#--output:sam bam bw mapping_picture
	#--need:(eval log); (eval pair); (eval work_dir(include plot_function.R));  map_stat.R
	fastq_dir=$1 
	name=$2
	reference_index=$3	
	out_dir=$4
	echo "the fastq is "$fastq_dir/$name".fq"
	echo "the reference" is $reference_index
	echo "the output_dir" is $out_dir
	hisat2 -p 10 -x $reference_index -1 $fastq_dir/$name"_R1.fq" -2 $fastq_dir/$name"_R2.fq" -S $out_dir/$name".sam" > $out_dir/$name".log" 2>&1 
	map_g=$(samtools view -h $out_dir/$name".sam" | grep -E 'NH:i' | awk '{print $1}' | sort | uniq | wc -l)
	echo $name"_map_genome "$map_g >> $log/$name"_record.txt"
	samtools view -bS $out_dir/$name".sam" > $out_dir/$name".bam"
	samtools sort -@ 8 $out_dir/$name".bam" -o $out_dir/$name"_sorted.bam"
	samtools index $out_dir/$name"_sorted.bam"
	rm -rf $out_dir/$name".bam"
	bamCoverage -b $out_dir/$name"_sorted.bam" -o $out_dir/$name".bw"
	mkdir $out_dir/tmp
	cat $out_dir/$name".log" | grep -E ') aligned concordantly' > $out_dir/tmp/tmp_mapstat.txt
	sed -i 's/[a-zA-Z]//g' $out_dir/tmp/tmp_mapstat.txt
	sed -i 's/[\(\)\%]//g' $out_dir/tmp/tmp_mapstat.txt
	sed -i 's/[ ][ ]*/ /g' $out_dir/tmp/tmp_mapstat.txt
	/usr/bin/Rscript map_stat.R -i $out_dir/tmp/tmp_mapstat.txt -w $out_dir/tmp -f $work_dir/"plot_function.R" -s $name > $out_dir/tmp/$name"_forstat.log" 2>&1 
	mv $out_dir/tmp/tmp_map.txt $out_dir/$name"_map.txt"
	mv $out_dir/tmp/tmp_mappie.pdf $out_dir/$name"_mappie.pdf"
	rm -rf $out_dir/tmp
}




















