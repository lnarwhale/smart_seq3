#!/bin/bash
function help_bulk(){
	echo "smartseq3_cut ｜目的:对smart数据进行质控｜位置:1.fastq所在文件夹 2.fastq的名字 3.质控后数据的输出文件夹"
	echo "---------------------------"
	echo "fastp_control | 目的:对数据进行fastp质控｜位置:1.fastq所在文件夹 2.fastq的名字 3.质控后数据的输出文件夹"
	echo "---------------------------"
	echo "fastqtobw | 目的:fastq比对,获取sam/bam/bw文件 | 位置:1.fastq所在文件夹 2.fastq的名字 3.map_index的路径 4.输出文件夹"
	echo "---------------------------"
	echo "rnaqc | 目的:对bam分析,从而获得降解/内外显子分布等基本信息 | 位置:1.bam所在文件夹 2.bam的名字 3.输出文件夹 4.注释gtf路径 5.注释bed路径"
	echo "---------------------------"
	echo "estimate_abundance | 目的:获取fpkm和rpkm | 位置:1.bam所在文件夹 2.bam的名字 3.注释gtf路径 4.输出文件夹"
	echo "---------------------------"
	echo "total_abundance | 目的:集合fpkm,并画箱型图 ｜位置:1.fpkm所在文件夹 2.输出文件夹"
	echo "---------------------------"
	echo "read_matrix | 目的:获取raw_coun和fpkm.csv ｜ 位置:1.工作路径 2.ballgown的路径" 
	echo "---------------------------"
	echo "trim_g | 目的:对fastq进行trim_galore得质控fastq ｜ 位置:1.fastq所在文件夹 2.fastq的名字 3.输出文件夹(不能与fastq所在文件夹一致) 4.trim的read长度阈值"
	echo "---------------------------"
	echo "de_rna | 目的:提取无法比对到指定index的read ｜ 位置:1.fastq所在文件夹 2.fastq名字 3.输出文件夹(不能与fastq所在文件夹一致); 4.指定index的路径"
	echo "---------------------------"
	echo "saturation | 目的:画饱和度 ｜ 位置:1.qualimeta.txt(rnaqc获取)所在文件夹 2.输出文件夹"
	echo "---------------------------"
	echo "TE_count | 目的:获取转座子的raw_count ｜ 位置:1.bam的名字 2.bam所在文件夹 3.输出文件夹 4.TE注释gtf的路径 5.注释gtf的路径"
	echo "---------------------------"
	echo "cager_a | 目的:对所有bam文件进行cager普通分析 | 位置:1.输出文件夹 2.物种 3.bam所在文件夹"
	echo "---------------------------"
	echo "tss_compare_single | 目的:对单端数据TSS进行edger差异分析 | 位置:1.compare(2c-1c|1c-mii) 2.out_dir 3.bam_dir 4.pwm_dir 5.species"	
	echo "---------------------------"
	echo "tss2loc | 目的:从cager_tss.csv文件获取基因symbol |位置:1.csv名字 2.csv_dir 3.gene_loc的gtf路径 4.out_dir"
	echo "---------------------------"
	echo "tss_TPM_compare | 目的:对数据进行TPM差异分析 ｜位置:1.compare组(2c-1c|1c-mii) 2.oue_dir 3.bam_dir 4.pwm_dir 5.species"
	echo "---------------------------"
	echo "cager_gene_tpm | 目的:经过promoter_c获取all_gene的tpm | 位置:1.cage工作文件夹 2.物种 3.bam_dir 4.gene_loc的gtf路径"
}
function smartseq3_cut(){
	#--function:cut the fastq of smart-seq3
	#--input: 1.fastq_dir; 2.name; 3.out_dir 
	#--output:fastq after cut
	#--need: precut.py; smartcut.py; (eval log)
	echo "the fastq is in "$1
	echo "the name is "$2
	echo "the output_dir is "$3
	mkdir $3/tmp
	raw=$(cat $1/$2".fq" | wc -l | awk '{print $1}')
	let raw=$raw/4
	/usr/bin/python precut.py -i $1 -o $3/tmp -n $2
	prec=$(cat $3/tmp/umi.fq | wc -l | awk '{print $1}')
	let prec=$prec/4
	/data1/shenluemou/biosoft/anaconda3/bin/seqkit rmdup --by-seq --ignore-case $3/tmp/umi.fq -D $3/$2"_umi.txt"> $3/tmp/deumi.fq
	all=$(wc -l $3/tmp/deumi.fq | awk '{print $1}')
	let all=$all/4	
	cat $3/$2"_umi.txt" | awk '{print $1}' > $3/$2"_umi_1.txt"
	chongfu=$(cat $3/$2"_umi_1.txt" | awk '{sum+=$1}END{print sum}')
	let dan=$all-$chongfu
	uniq -c $3/$2"_umi_1.txt" > $3/$2"_umi_2.txt"
	echo " "$dan" 1" >> $3/$2"_umi_2.txt"
	sed -i 's/ \+/ /g' $3/$2"_umi_2.txt"
	sed -i 's/ /,/g' $3/$2"_umi_2.txt"
	sed -i  's/^.//g' $3/$2"_umi_2.txt"
	mv $3/$2"_umi_2.txt" $log/$2"_umi_sum.txt"
	nodu=$(cat $3/tmp/deumi.fq | wc -l | awk '{print $1}')
	let nodu=$nodu/4
	/usr/bin/python smartcut.py -i $3/tmp -o $3 -n $2
	rm -rf $3/tmp
	touch $log/$2"_record.txt"
	echo $2"_raw "$raw >> $log/$2"_record.txt"
	echo $2"_prec "$prec >> $log/$2"_record.txt"
	echo $2"_nodu "$nodu >> $log/$2"_record.txt"
}
#----------------------------------------------------------------------
function fastp_control(){
	#--function:fastp quality control
	#--input:1.fastq_dir 2.name 3.out_dir
	#--need:fastp
	echo "the fastp of $1/$2 begin"
	/data1/shenluemou/biosoft/fastp \
	-i $1/$2".fq" \
	-o $3/$2".fq" \
	--html $3/$2".html" \
	--trim_poly_g --poly_g_min_len 5 \
	--trim_poly_x --poly_x_min_len 5 \
	--cut_front --cut_tail --cut_window_size 4 \
	--qualified_quality_phred 15 \
	--low_complexity_filter \
        --complexity_threshold 30 \
        --length_required 30 \
        --thread 4
}
#-------------------------------------------------------------------
function fastqtobw(){
	#--function:to get bw and mapping_picture from fastq
	#--input:1.fastq_dir; 2.name; 3.reference_dir&name(index); 4.out_dir 
	#--output:sam bam bw mapping_picture
	#--need: (eval log); (eval pair); (eval work_dir(include plot_function.R));  map_stat.R
	echo "the fastq is in "$1
	echo "the name is "$2
	echo "the reference is "$3
	echo "the output_dir is "$4
	hisat2 -p 10 -x $3 -U $1/$2".fq" -S $4/$2".sam" > $4/$2".log" 2>&1
	map_g=$(samtools view -h $4/$2".sam" | grep -E 'NH:i' | awk '{print $1}' | sort | uniq | wc -l)
	echo $2"_map_genome "$map_g >> $log/$2"_record.txt" 
	samtools view -bS $4/$2".sam" > $4/$2".bam"
	samtools sort -@ 8 $4/$2".bam" -o $4/$2"_sorted.bam" 
	samtools index $4/$2"_sorted.bam"
	rm -rf $4/$2".bam"
	bamCoverage -b $4/$2"_sorted.bam" -o $4/$2".bw"
	mkdir $4/tmp
	if [ $pair = "double" ];then
		cat $4/$2".log" | grep -E ') aligned concordantly' > $4/tmp/tmp_mapstat.txt
	fi
	if [ $pair = "single" ];then
		cat $4/$2".log" | grep -E 'aligned' > $4/tmp/tmp_mapstat.txt
	fi
	sed -i 's/[a-zA-Z]//g' $4/tmp/tmp_mapstat.txt
	sed -i 's/[\(\)\%]//g' $4/tmp/tmp_mapstat.txt
	sed -i 's/[ ][ ]*/ /g' $4/tmp/tmp_mapstat.txt
	/usr/bin/Rscript map_stat.R -i $4/tmp/tmp_mapstat.txt -w $4/tmp -f $work_dir/"plot_function.R" -s $2 > $4/tmp/$2"_forstat.log" 2>&1
	mv $4/tmp/tmp_map.txt $4/$2"_map.txt"
	mv $4/tmp/tmp_mappie.pdf $4/$2"_mappie.pdf"
	rm -rf $4/tmp
}
#-----------------------------------------------------------------
function rnaqc(){
	#--function:to do the rna_qc for bam
	#--input:1.bam_dir; 2.name; 3.out_dir; 4.reference_dir&name(gtf);	5.reference_dir&name(bed)
	#--output:picture of distribution and picture of deletation
	#--need: (eval log); (eval work_dir(include plot_function.R)); qualimap; quali.R; read_distribution.py; rseqc_stat.R
	mkdir $3/log
	mkdir $3/tmp
	mkdir $3/read_distribution/
	echo "the bam's name is "$2
	echo "the bam is in "$1
	echo "the output_dir is "$3
	echo "gtf is "$4
	echo "bed is "$5
	./qualimap/qualimap --java-mem-size=12G rnaseq -bam $1/$2"_sorted.bam" -gtf $4 -outdir $3/tmp/ -pe -outformat PDF:HTML > $3/$2".log" 2>&1
	cat $3/tmp/rnaseq_qc_results.txt | grep -E 'exonic|intronic|intergenic' > $3/tmp/tmp_readdis.txt
	sed -i 's/[\=\(\)\%\,]//g' $3/tmp/tmp_readdis.txt
	sed -i 's/[ ][ ]*/ /g' $3/tmp/tmp_readdis.txt
	sed -i 's/[\#]//g' $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(total).txt"
	sed -i 's/[\#]//g' $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(high).txt"
	sed -i 's/[\#]//g' $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(low).txt"	
	/usr/bin/Rscript quali.R -d $3/tmp/tmp_readdis.txt -t $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(total).txt" -h $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(high).txt" -l $3/tmp/raw_data_qualimapReport/"coverage_profile_along_genes_(low).txt" -w $3/tmp -f ${work_dir}/"plot_function.R" -s $2 > $3/tmp/$2"_quali.log" 2>&1
	mv $3/tmp/tmp_quali.txt  $3/read_distribution/$2"_quali.txt"
	mv $3/tmp/tmp_readdis_pie.pdf $3/read_distribution/$2"_readdis_pie.pdf"
	mv $3/tmp/tmp_readdis.pdf $3/read_distribution/$2"_readdis_line.pdf"
	mkdir $3/read_distribution/tmp/
	/home/shenlm/.local/bin/read_distribution.py -i $1/$2"_sorted.bam" -r $5 > $3/read_distribution/tmp/tmp_tag.txt
	cp $3/read_distribution/tmp/tmp_tag.txt $3/read_distribution/$2"_taq.txt"
	cat $3/read_distribution/$2"_taq.txt" | grep -E 'Group|Exons|Introns|TSS|TES' > $3/read_distribution/tmp/tag.txt
	sed -i 's/[ ][ ]*/ /g' $3/read_distribution/tmp/tag.txt
	sed -i  "s/[\']//g" $3/read_distribution/tmp/tag.txt
	/usr/bin/Rscript rseqc_stat.R -i $3/read_distribution/tmp/tag.txt -w $3/read_distribution/tmp/ -f ${work_dir}/"plot_function.R" -s $2 > $3/read_distribution/tmp/$2"reddis.log" 2>&1
	mv $3/read_distribution/tmp/tmp_rseqc.txt $3/read_distribution/$2"_rseqc.txt"
	mv $3/read_distribution/tmp/tmp_readdis.pdf $3/read_distribution/$2"_readdis_bar.pdf"
	rm -rf $3/read_distribution/tmp/
}
#-----------------------------------------------------------
function estimate_abundance(){
	#--function: using stringtie to calculate the rpkm or fpkm
	#--input:1.bam_dir; 2.name; 3.reference_dir&name(gtf) 4.output_dir
	#--output: fpkm/rpkm.txt
	#--need: (eval log); (eval work_dir(include plot_function.R))
	echo "the bam is"$1"/"$2"_sorted.bam"
	echo "the gtf is "$3
	echo "the output_dir is "$4
	mkdir $4/ballgown/
	mkdir $4/ballgown/$2
	stringtie -e -B -A $4/$2 -p 8 -G $3 -o $4/ballgown/$2/$2"_str_quant.gtf" $1/$2"_sorted.bam"
	awk '{print $1,$8}' $4/$2 > $4/$2"_fpkm_tmp.txt"
	sed -i '1d' $4/$2"_fpkm_tmp.txt"
	awk '$0=$0" '$2'"' $4/$2"_fpkm_tmp.txt" > $4/$2"_fpkm.txt"
	rm -rf $4/$2"_fpkm_tmp.txt"
}
#----------------------------------------------------------
function total_abundance(){
	#--function: merge fpkm and box-plot
	#--input: 1.fpkm_dir 2.working_dir
	#--output: fpkm_all.csv ; box-plot
	echo "-------------------------------------"
	echo "the fpkm is in "$1
	echo "the working_dir/output-dir is in "$2
	cat $1/*.txt > $2/total_expression.txt
	/usr/bin/Rscript fpkm_dis.R -i $2/total_expression.txt -w $2 -f ${work_dir}/plot_function.R > $2/fpkm_dis.log 2>&1
	echo "the total_abundance finished"
	echo "-------------------------------------"
}
#----------------------------------------------------------
function read_matrix(){
	#--funtion:generate box_plot and read_matrix
	#--input: 1.working_dir 2.ballgown_dir 
	#--output: gene_count.csv ; fpkm.csv
	#--need : eval work_dir ; split_countmatrix.R
	echo "the working_dir is "$1
	echo "the ballgown_dir is "$2
	cd $1
	cp -r $2/ballgown $1
	/data1/shenluemou/con/prepDE.py -i ./ballgown > ./read_matrix.log 2>&1
	cd ${work_dir}
	/usr/bin/Rscript split_countmatrix.R -i $1/gene_count_matrix.csv -w $1 >> ./read_matrix.log 2>&1
	echo "the ballgown has finished"
	echo "---------------------------------------------------------------------------------"
}
#----------------------------------------------------------
function trim_g(){
	#--function:trim_galore
	#--input: 1.fastq_dir; 2.name; 3.output_dir(must different from fastq_dir) 4.trim_length
	#--output:fastq after trim
	#--need:(eval pair)
	echo "the fastq is "$1"/"$2".fq"
	echo "the output is "$3
	echo "the pair is "${pair}
	if [ $pair = "single" ];then
		/data1/shenluemou/biosoft/TrimGalore-0.6.6/trim_galore -q 25 --phred33 --stringency 3 --length $4 -e 0.1 --no_report_file $1/$2".fq" -o $3/tmp/ > $3/$2".log" 2>&1
		mv $3/tmp/$2"_trimmed.fq" $3/$2".fq"
		rm -rf $3/tmp/
		tri=$(cat $3/$2".fq"|wc -l|awk '{print $1}')
		let tri=$tri/4
		echo $2"_trim "$tri >> $log/$2"_record.txt"
    fi
    if [ $pair = "double" ];then
        /data1/shenluemou/biosoft/TrimGalore-0.6.6/trim_galore -q 25 --phred33 --stringency 3 --length $4 -e 0.1 --no_report_file --paired $1/$2"_R1.fq" $1/$2"_R2.fq" -o $3/tmp/ > $3/$2".log" 2>&1
		mv $3/tmp/$2"_R1_trimmed.fq" $3/$2"_R1.fq"
		mv $3/tmp/$2"_R2_trimmed.fq" $3/$2"_R2.fq"
		rm -rf $3/tmp/
    fi
}
#---------------------------------------------------------
function de_rna(){
	#--function:delete the rna 
	#--input:1.fastq_dir; 2.name; 3.out_dir(different from input); 4.rna_index_dir&name(index)
	#--output:plot of rna%; fastq afer delete rna
	#--need:(eval work_dir);map_stat.R
	hisat2 -p 10 -x $4 -U $1/$2".fq" -S $3/$2"_rna.sam" > $3/$2".log" 2>&1
	cat $3/$2"_rna.sam" | grep -v "NH:i" > $3/$2"_norna.sam" 
	samtools fastq $3/$2"_norna.sam" > $3/$2".fq"
	de_rrna=$(wc -l $3/$2".fq"| awk '{print $1}')
	let de_rrna=$de_rrna/4
	echo  $2"_map_rrna "$de_rrna >> $log/$2"_record.txt"
	rm -rf $3/$2"_norna.sam" $3/$2"_rna.sam"
	mkdir $3/tmp
	cat $3/$2".log" | grep -E 'aligned' > $3/tmp/tmp_mapstat.txt
	sed -i 's/[a-zA-Z]//g' $3/tmp/tmp_mapstat.txt
	sed -i 's/[\(\)\%]//g' $3/tmp/tmp_mapstat.txt
	sed -i 's/[ ][ ]*/ /g' $3/tmp/tmp_mapstat.txt
	/usr/bin/Rscript map_stat.R -i $3/tmp/tmp_mapstat.txt -w $3/tmp -f $work_dir/"plot_function.R" -s $2 > $3/tmp/$2"_forstat.log" 2>&1
	mv $3/tmp/tmp_map.txt $3/$2"_map.txt"
	mv $3/tmp/tmp_mappie.pdf $3/$2"_mappie.pdf"
}
#----------------------------------------------------------------
function extract_high(){
	#--function: extract the high-expression gene from sample according control
	#--input: 1.rawcount_dir&name 2.outdir 3.control_name 4.sample_name
	#--output: high_expression gene
	/usr/bin/Rscript extrac_high_expression.R -t $1 -w $2 -c $3 -s $4
}
#---------------------------------------------------------------
function saturation(){
	#--function:saturation plot
	#--input: 1.qualimeta.txt(absolute path and include "qualimeta.txt") 2.out_dir
	#--output:saturation plot
	./qualimap/qualimap counts -d $1 -outdir $2 >> $2/"read_matrix.log"
}
#---------------------------------------------------------------
function TE_count(){
	#--function:TE count
	#--input: 1.name 2.bam_dir 3.outdir 4.TE_reference_dir&name 5.gtf_dir&name
	#--outdir:heatmap of TE count;txt of TE count
	TEcount --mode multi \
	-b $2/$1"_sorted.bam" \
	--GTF $5 \
	--TE $4 \
	--project $3/$1 \
	--sortByPos
	cat $3/$1".cntTable" | grep -v "ENS" | sed '1d' | sed 's/\t/,/g' | sed 's/:/,/g' > $3/$1"_te.csv"
}
#---------------------------------------------------------------
function cager_a(){
	#--input:1.out_dir 2.species 3.bam_dir(*_sorted.bam)
	#--output:1.promoter_picture 2.width_picture 
	Rscript cager_all.R -w $1 -s $2 -b $3
}
#---------------------------------------------------------------
function tss_compare_single(){
	#--input:1.compare组 2.out_dir 3.bam_dir 4.pwm_dir 5.species
	#--ouput:1.差异tss位点 2.upset 3.tata-tct饼图
	compare=$1
	species=$5
	compare_tion="1 2"
	bam_dir=$3
	all_dir=$2
	pwm_dir=$4
	compare=${compare//|/ }
	for condition in ${compare};do
        	mkdir $all_dir/$condition
        	compare_dir=$all_dir/$condition
        	case="${condition%-*}"
        	control="${condition#*-}"
        	cp $bam_dir/$case"_sorted.bam" $compare_dir/"x2_"$case"_sorted.bam"
        	cp $bam_dir/$control"_sorted.bam" $compare_dir/"x1_"$control"_sorted.bam"
        	for sec in ${compare_tion};do
                	mkdir $compare_dir/result
                	Rscript cager_edger_set.R -s $species -o $compare_dir/"result" -b $compare_dir -p $pwm_dir -m $sec
        	done
	done
}
#-------------------------------------------------------------
function tss2loc(){
	#--input:1.name 2.csv_dir 3.gene_loc_gtf_index 4.out_dir
	#--output:gene_loc.csv
	name=$1
	indir=$2
	gene_loc_gtf=$3
	outdir=$4
	cat $indir/$name".csv" | cut -d , -f 2,3,4 | sed 's/"//g' | sed 's/,/\t/g' | sed '1d' > $outdir/$name"_x1_loc.csv"
	bedtools intersect -a $outdir/$name"_x1_loc.csv" -b $gene_loc_gtf -wa -wb | bedtools groupby -i - -g 1-4 -c 7 -o collapse > $outdir/$name"_x2_locgene.csv"
	cat $outdir/$name"_x2_locgene.csv" | sed 's/\t/,/g' | cut -d , -f 5 | sort | uniq > $outdir/$name"_x3_gene.csv"
}
#-------------------------------------------------------------
function tss_TPM_compare(){
	#--input:1.compare组 2.oue_dir 3.bam_dir 4.pwm_dir 5.species
	#--outpit:1.差异tss位点 2.upset 3.tata-tct饼
	compare=$1
	species=$5
	bam_dir=$3
	all_dir=$2
	pwm_dir=$4
	compare=${compare//|/ }
	for condition in ${compare};do
		mkdir $all_dir/$condition
		case="${condition%-*}"
		control="${condition#*-}"
		compare_dir=$all_dir/$condition
		cp $bam_dir/$case"_sorted.bam" $compare_dir/"x2_"$case"_sorted.bam"
		cp $bam_dir/$control"_sorted.bam" $compare_dir/"x1_"$control"_sorted.bam"
		mkdir $compare_dir/result
		Rscript cager_TPM_diff.R -s $species -o $compare_dir/"result" -b $compare_dir -p $pwm_dir  
	done
}
#-------------------------------------------------------------
function cager_gene_tpm(){
	#--input:
	#--output:
	work_dir=$1
	species=$2
	bam_dir=$3
	gene_loc=$4
	Rscript cager_gene.R -w $work_dir -s $species -b $bam_dir
	cat $work_dir/names.csv | cut -d , -f 2 | sed 's/"//g' | sed '1d'| tr '\n' ' ' > $work_dir/names_all.csv
	group=$(cat $work_dir/names_all.csv)
	mkdir $work_dir/tss
	for name in ${group};do
		tss2loc $name"_cager_tss" $work_dir $gene_loc $work_dir/tss
	done
}












