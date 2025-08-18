#!/bin/sh
#SBATCH --job-name=job
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4GB

module load python-anaconda3
source activate rarjun

DATADIR="../fastqs/BulkRNA/*"
FASTQ_SCREEN_GENOMES="/lsi/home/rarjun/FastQ_Screen_Genomes/fastq_screen.conf" 


for f in $DATADIR/*.fastq
 	do
 		pigz $f
	done;

for d in $DATADIR
do
 	folder_name="${d%/}"
    
    if [ -z "$(ls -A "$folder_name")" ]; then
		echo "Skipping empty folder: $folder_name"
		continue
    fi
	experiment=$(basename $d)

 	mkdir ../QC/${experiment}
    chmod 777 ../QC/${experiment}
	
    mkdir ../processed/${experiment}
    chmod 777 ../processed/${experiment}

	mkdir ../peaks/${experiment}
	chmod 777 ../peaks/${experiment}
 
	echo $experiment
	query=$(find ${d}/ -type f -name "*1_001.fastq.gz")
	files=($query)
	FASTQ1=${files[0]}
	query=$(find ${d}/ -type f -name "*2_001.fastq.gz")
	files=($query)
	FASTQ2=${files[0]}

    
	fastqc -o ../QC/${experiment} $FASTQ1 $FASTQ2 -t 16
	fastq_screen --threads 16 --conf $FASTQ_SCREEN_GENOMES --outdir ../QC/${experiment} $FASTQ1 $FASTQ2 
    
	BAM="../processed/"${experiment}"/"${experiment}"_sorted.bam"

	if [ ! -e "$BAM" ]; then
		echo "Bam file doesnt exist, processing"
		cutadapt -j 16 \
			-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
			-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
			-o ../processed/${experiment}/R1.trimmed.fastq \
			-p ../processed/${experiment}/R2.trimmed.fastq \
			$FASTQ1 $FASTQ2 --json="../QC/${experiment}/${experiment}_cutadapt_report.json"

		hisat2 -q -t -x genome_tran \
			-1 ../processed/${experiment}/R1.trimmed.fastq \
			-2 ../processed/${experiment}/R2.trimmed.fastq \
			-S ../processed/${experiment}/aligned.sam \
			--no-mixed --no-discordant --quiet --no-unal \
			--summary-file "../QC/${experiment}/${experiment}_hisat_report.txt" -p 16

		samtools view -bS ../processed/${experiment}/aligned.sam | samtools sort - -o $BAM -@ 16
		samtools index $BAM 
	else
		echo "Bam file exists. Skipping"
	fi
done


featureCounts -a "dm6.refGene.gtf" \
	-o ../analysis_data/bulkRNA/RNA_featureCounts.txt \
	-F 'GTF' \
	-p \
	-T 64 \
	../processed/*RNA*/*.bam