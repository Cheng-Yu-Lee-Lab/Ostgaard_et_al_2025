#!/bin/sh
#SBATCH --job-name=job
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4GB

function run_gopeaks() {
file=$1
IgG=$2
sample=$3
sbatch  --job-name=${3} --output=${3}.out --error=${3}.err -c 6 --mem-per-cpu=4GB <<END
#!/bin/bash
gopeaks -b $1 -c $2 -o ../peaks_CR/gopeaks/broad/$3 --chromsize dm6.genome --broad
gopeaks -b $1 -c $2 -o ../peaks_CR/gopeaks/narrow/$3 --chromsize dm6.genome  
END
}

mkdir ../peaks_CR/gopeaks
mkdir ../peaks_CR/gopeaks/broad
mkdir ../peaks_CR/gopeaks/narrow
IgG_file="../processed_CR/IgG/IgG_sorted_dedup.bam"
for f in ../processed_CR/H3K*/*.bam
do
    sample=$(basename "$f" | cut -d_ -f1-3)
    firstCharacter=${sample:0:1}
    H="H"
    K="K"
    
    if [ "$firstCharacter" == "$H" ] || [ "$firstCharacter" == "$K" ]; then
		echo $sample
		run_gopeaks $f $IgG_file $sample
    fi 
done

mkdir ../peaks_CR/gopeaks/merged
mkdir ../peaks_CR/gopeaks/blacklisted
mkdir ../peaks_CR/gopeaks/SAF

cat ../peaks_CR/gopeaks/broad/H3K27ac_*_peaks.bed \
  | sort -k1,1 -k2,2n | bedtools merge -i stdin > ../peaks_CR/gopeaks/merged/H3K27ac_broad.bed
  
cat ../peaks_CR/gopeaks/broad/H3K27me3_*_peaks.bed \
  | sort -k1,1 -k2,2n | bedtools merge -i stdin > ../peaks_CR/gopeaks/merged/H3K27me3_broad.bed
 
cat ../peaks_CR/gopeaks/narrow/H3K4me1_*_peaks.bed \
  | sort -k1,1 -k2,2n | bedtools merge -i stdin > ../peaks_CR/gopeaks/merged/H3K4me1_narrow.bed

cat ../peaks_CR/gopeaks/narrow/H3K4me3_*_peaks.bed \
   | sort -k1,1 -k2,2n | bedtools merge -i stdin > ../peaks_CR/gopeaks/merged/H3K4me3_narrow.bed

for f in ../peaks_CR/gopeaks/merged/*.bed
do
     temp=${f#*/}
     sample=$(basename "$f" .bed)
     echo $sample
     bedtools intersect -a $f -b dm6.blacklist.noCHR.bed -v > ../peaks_CR/gopeaks/blacklisted/${sample}.blacklisted.bed
     awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' ../peaks_CR/gopeaks/blacklisted/${sample}.blacklisted.bed > ../peaks_CR/gopeaks/SAF/${sample}.saf
done

mkdir ../CR_counts
featureCounts -a ../peaks_CR/gopeaks/SAF/H3K27ac_broad.saf \
	-o ../CR_counts/H3K27ac.txt \
	-F 'SAF' \
	-p \
	--countReadPairs \
	../processed_CR/H3K27ac*/*.bam
	
featureCounts -a ../peaks_CR/gopeaks/SAF/H3K27me3_broad.saf \
	-o ../CR_counts/H3K27me3.txt \
	-F 'SAF' \
	-p \
	--countReadPairs \
	../processed_CR/H3K27me3*/*.bam
	
featureCounts -a ../peaks_CR/gopeaks/SAF/H3K4me1_narrow.saf \
	-o ../CR_counts/H3K4me1.txt \
	-F 'SAF' \
	-p \
	--countReadPairs \
	../processed_CR/H3K4me1*/*.bam

featureCounts -a ../peaks_CR/gopeaks/SAF/H3K4me3_narrow.saf \
	-o ../CR_counts/H3K4me3.txt \
	-F 'SAF' \
	-p \
	--countReadPairs \
	../processed_CR/H3K4me3*/*.bam
	
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' \
  ../analysis_data/TF_peaks/ase.bed > ../analysis_data/TF_peaks/ase.saf

awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' \
  ../analysis_data/TF_peaks/btd.bed > ../analysis_data/TF_peaks/btd.saf

awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' \
  ../analysis_data/TF_peaks/fruC-myc.bed > ../analysis_data/TF_peaks/fruC-myc.saf

awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' \
  ../analysis_data/TF_peaks/fruCOM_aPKC.bed > ../analysis_data/TF_peaks/fruCOM_aPKC.saf

awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' \
  ../analysis_data/TF_peaks/zld.bed > ../analysis_data/TF_peaks/zld.saf

awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print NR, $1, $2+1, $3, "."}' \
  ../analysis_data/TF_peaks/zld_aPKC.bed > ../analysis_data/TF_peaks/zld_aPKC.saf

featureCounts -a ../analysis_data/TF_peaks/ase.saf \
  -o ../CR_counts/ase.txt \
  -F 'SAF' \
  -p \
  --countReadPairs \
  ../processed_CR/*ase*/*.bam

featureCounts -a ../analysis_data/TF_peaks/btd.saf \
  -o ../CR_counts/btd.txt \
  -F 'SAF' \
  -p \
  --countReadPairs \
  ../processed_CR/*btd*/*.bam

featureCounts -a ../analysis_data/TF_peaks/fruC-myc.saf \
  -o ../CR_counts/fruC.txt \
  -F 'SAF' \
  -p \
  --countReadPairs \
  ../processed_CR/*FruC_0hr*/*.bam

featureCounts -a ../analysis_data/TF_peaks/fruCOM_aPKC.saf \
  -o ../CR_counts/fruCOM_aPKC.txt \
  -F 'SAF' \
  -p \
  --countReadPairs \
  ../processed_CR/*FruCOM_apkc*/*.bam

featureCounts -a ../analysis_data/TF_peaks/zld.saf \
  -o ../CR_counts/zld.txt \
  -F 'SAF' \
  -p \
  --countReadPairs \
  ../processed_CR/*Zld_0hr*/*.bam

featureCounts -a ../analysis_data/TF_peaks/zld_aPKC.saf \
  -o ../CR_counts/zld_aPKC.txt \
  -F 'SAF' \
  -p \
  --countReadPairs \
  ../processed_CR/*Zld_apkc*/*.bam