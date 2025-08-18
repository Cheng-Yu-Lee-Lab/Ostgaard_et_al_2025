#!/bin/sh
#SBATCH --job-name=job
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4GB

module load python-anaconda3
source activate rarjun

DATADIR="../fastqs/CUT_RUN/*"
FASTQ_SCREEN_GENOMES="/lsi/home/rarjun/FastQ_Screen_Genomes/fastq_screen.conf" 


function submit_job() {
experiment=$1
d=$2
fastq=$3
FASTQ1=( "$d"/"$fastq"R1_001.fastq.gz )
FASTQ2=( "$d"/"$fastq"R2_001.fastq.gz )
FASTQ1_gz=( "$d"/"$fastq"R1_001.fastq )
FASTQ2_gz=( "$d"/"$fastq"R2_001.fastq )

BAM="../processed_CR/${experiment}/${experiment}_sorted_dedup.bam"
#--mail-user=rarjun@umich.edu --mail-type=ALL
sbatch  --job-name=${1} --output=${1}.out --error=${1}.err -c 6 --mem-per-cpu=8GB <<END
#!/bin/bash
echo $d

pigz $FASTQ1_gz
pigz $FASTQ2_gz

echo $FASTQ1
echo $FASTQ2

mkdir ../QC_CR/${experiment}
chmod 777 ../QC_CR/${experiment}

mkdir ../processed_CR/${experiment}
chmod 777 ../processed_CR/${experiment}

mkdir ../peaks_CR/${experiment}
chmod 777 ../peaks_CR/${experiment}

echo $FASTQ1
echo $FASTQ2
whereis bowtie2

fastqc -o ../QC_CR/${experiment} $FASTQ1 $FASTQ2 -t 32
fastq_screen --threads 32 --conf $FASTQ_SCREEN_GENOMES --outdir ../QC_CR/${experiment} $FASTQ1 $FASTQ2 

BAM="../processed_CR/"${experiment}"/"${experiment}"_sorted_dedup.bam"
echo $BAM

if [ ! -e "$BAM" ]; then
echo "Bam file doesnt exist, processing"
NGmerge -a -1 $FASTQ1 -2 $FASTQ2 -o ../processed_CR/${experiment}/NGMerge -v -n 6

bowtie2 -x BDGP6 \
	 -1 ../processed_CR/${experiment}/NGMerge_1.fastq.gz -2 ../processed_CR/${experiment}/NGMerge_2.fastq.gz \
	 --very-sensitive -I 10 -X 2000 \
	 --no-unal -p 6 \
	 -S ../processed_CR/${experiment}/${experiment}_aligned.sam \
	 > ../QC_CR/${experiment}/${experiment}.bowtie2_error.log 2> ../QC_CR/${experiment}/${experiment}.bowtie2_aligned.log

samblaster -i ../processed_CR/${experiment}/${experiment}_aligned.sam --removeDups | samtools view -bS - | samtools sort - -o $BAM -@ 32
else
# Code to run when the file exists
echo "Bam file exists. Skipping"
fi

samtools index $BAM 
rm ../processed_CR/${experiment}/${experiment}_aligned.sam

bamCoverage -b $BAM -o ../bigwigs_CR/${experiment}.small.bw --binSize 5 -p max --normalizeUsing RPKM --maxFragmentLength 120
bamCoverage -b $BAM -o ../bigwigs_CR/${experiment}.bw --binSize 5 -p max --normalizeUsing RPKM

macs2 callpeak -t $BAM -g 142573017 --outdir ../peaks_CR/${experiment} -n ${experiment}

END
}

for d in $DATADIR
do
   experiment=$(basename $d)
	first_file=$(find "$d" -maxdepth 1 -type f -printf '%f\n' | sort | head -n1)
	fastq="${first_file%%R1_001.*}"
	echo $fastq
   submit_job $experiment $d $fastq
done

