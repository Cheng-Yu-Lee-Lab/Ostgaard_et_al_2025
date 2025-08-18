#!/bin/sh
#SBATCH --job-name=job
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4GB

module load python-anaconda3
source activate rarjun

DATADIR="../fastqs/ATAC/*"
FASTQ_SCREEN_GENOMES="/lsi/home/rarjun/FastQ_Screen_Genomes/fastq_screen.conf" 

mkdir ../QC_ATAC/atacqv_results/
chmod 777 ../QC_ATAC/atacqv_results/

for d in $DATADIR
do
	experiment=$(basename $d)

	mkdir ../QC_ATAC/${experiment}
    chmod 777 ../QC_ATAC/${experiment}

    mkdir ../processed_ATAC/${experiment}
    chmod 777 ../processed_ATAC/${experiment}

    mkdir ../peaks_ATAC/${experiment}
    chmod 777 ../peaks_ATAC/${experiment}
 
    echo $experiment
    query=$(find ${d}/ -type f -name "*1_001.fastq")
    files=($query)
    FASTQ1=${files[0]}
    query=$(find ${d}/ -type f -name "*2_001.fastq")
    files=($query)
    FASTQ2=${files[0]}

    echo $FASTQ1
    echo $FASTQ2
    whereis bowtie2
    
    fastqc -o ../QC_ATAC/${experiment} $FASTQ1 $FASTQ2 -t 32
    fastq_screen --threads 32 --conf $FASTQ_SCREEN_GENOMES --outdir ../QC_ATAC/${experiment} $FASTQ1 $FASTQ2 
    
    BAM="../processed_ATAC/"${experiment}"/"${experiment}"_sorted_dedup.bam"
    echo $BAM

    if [ ! -e "$BAM" ]; then
    	echo "Bam file doesnt exist, processing"
    	NGmerge  -a  -1 $FASTQ1 -2 $FASTQ2 -o ../processed_ATAC/${experiment}/NGMerge -v -n 6

	    bowtie2 -x BDGP6 \
	         -1 ../processed_ATAC/${experiment}/NGMerge_1.fastq -2 ../processed_ATAC/${experiment}/NGMerge_2.fastq \
	         --very-sensitive -I 10 -X 2000 \
	         --no-unal -p 6 \
	         -S ../processed_ATAC/${experiment}/${experiment}_aligned.sam \
	         > ../QC_ATAC/${experiment}/${experiment}.bowtie2_error.log 2> ../QC_ATAC/${experiment}/${experiment}.bowtie2_aligned.log
	    samblaster -i ../processed_ATAC/${experiment}/${experiment}_aligned.sam --removeDups | samtools view -bS - | samtools sort - -o $BAM -@ 32
	else
    	echo "Bam file exists. Skipping"
	fi
    
    samtools index $BAM 
    bamCoverage -b $BAM -o ../bigwigs_ATAC/${experiment}.small.bw --binSize 5 -p max --normalizeUsing RPKM --maxFragmentLength 120
	macs2 callpeak -t $BAM -g 142573017 --outdir ../peaks_ATAC/${experiment} -n ${experiment}
	
	ataqv --peak-file ../peaks_ATAC/${experiment}_peaks.narrowPeak \
	  --name ${experiment} \
	  --metrics-file ../QC_ATAC/atacqv_results/${experiment}.json.gz \
	  --excluded-region-file dm6.blacklist.noCHR.bed \
	  --tss-file dm6.tss.bed \
	  fly \
	  $BAM > ../QC_ATAC/atacqv_results/${experiment}.ataqv.out
 done

mkarv ../QC_ATAC/atac_comp ../QC_ATAC/atacqv_results/* --force

cat ../peaks_ATAC/*/*.narrowPeak \
  | sort -k1,1 -k2,2n | bedtools merge -i stdin | \
  bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v \
  > ../peaks_ATAC/merged_ATAC.blacklisted.bed
  
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' ../peaks_ATAC/merged_ATAC.blacklisted.bed > ../peaks_ATAC/merged_ATAC.blacklisted.saf

mkdir ../ATAC_counts

featureCounts -a ../peaks_ATAC/merged_ATAC.blacklisted.saf \
	-o ../ATAC_counts/ATAC_counts.txt \
	-F 'SAF' \
	-p \
	--countReadPairs \
	../processed_ATAC/*/*.bam