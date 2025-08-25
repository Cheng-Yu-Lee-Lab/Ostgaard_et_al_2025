#!/bin/sh
#SBATCH --job-name=job
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4GB

mkdir ../bigwigs_ATAC/TMM_norm

bamCoverage --bam ../processed_ATAC/0hr_rep1/*.bam -o ../bigwigs_ATAC/TMM_norm/0hr_rep1.small.bw --binSize 5 --scaleFactor 0.06587381 --numberOfProcessors 16 --maxFragmentLength 120
bamCoverage --bam ../processed_ATAC/0hr_rep2/*.bam -o ../bigwigs_ATAC/TMM_norm/0hr_rep2.small.bw --binSize 5 --scaleFactor 0.06315934 --numberOfProcessors 16 --maxFragmentLength 120
bamCoverage --bam ../processed_ATAC/0hr_rep3/*.bam -o ../bigwigs_ATAC/TMM_norm/0hr_rep3.small.bw --binSize 5 --scaleFactor 0.07063151 --numberOfProcessors 16 --maxFragmentLength 120
bamCoverage --bam ../processed_ATAC/24hr_rep1/*.bam -o ../bigwigs_ATAC/TMM_norm/24hr_rep1.small.bw --binSize 5 --scaleFactor 0.07331903 --numberOfProcessors 16 --maxFragmentLength 120
bamCoverage --bam ../processed_ATAC/24hr_rep2/*.bam -o ../bigwigs_ATAC/TMM_norm/24hr_rep2.small.bw --binSize 5 --scaleFactor 0.07737281 --numberOfProcessors 16 --maxFragmentLength 120
bamCoverage --bam ../processed_ATAC/24hr_rep3/*.bam -o ../bigwigs_ATAC/TMM_norm/24hr_rep3.small.bw --binSize 5 --scaleFactor 0.08085007 --numberOfProcessors 16 --maxFragmentLength 120
bamCoverage --bam ../processed_ATAC/apkc_rep1/*.bam -o ../bigwigs_ATAC/TMM_norm/apkc_rep1.small.bw --binSize 5 --scaleFactor 0.04306116 --numberOfProcessors 16 --maxFragmentLength 120
bamCoverage --bam ../processed_ATAC/apkc_rep2/*.bam -o ../bigwigs_ATAC/TMM_norm/apkc_rep2.small.bw --binSize 5 --scaleFactor 0.04770262 --numberOfProcessors 16 --maxFragmentLength 120
bamCoverage --bam ../processed_ATAC/apkc_rep3/*.bam -o ../bigwigs_ATAC/TMM_norm/apkc_rep3.small.bw --binSize 5 --scaleFactor 0.04984023 --numberOfProcessors 16 --maxFragmentLength 120

mkdir ../bigwigs_ATAC/TMM_merged
INDIR="../bigwigs_ATAC/TMM_norm/"
OUTDIR="../bigwigs_ATAC/TMM_merged"

for rep1 in $INDIR/*_rep1.small.bw; do
  base=$(basename "$rep1")                 
  core="${base%_rep1.small.bw}"              
  rep2="$INDIR/${core}_rep2.small.bw"         
  rep3="$INDIR/${core}_rep3.small.bw"         
  out="$OUTDIR/${core}.bg"                

  if [[ -f "$rep2" ]]; then
    echo "Merging: $base  +  $(basename "$rep2")  ->  $(basename "$out")"
	wiggletools write_bg "$out" mean \
                 "$rep1" \
                 "$rep2" \
				 "$rep3"				 
  else
    echo "WARNING: Missing pair for ${core}: ${rep2} not found. Skipping." >&2
  fi
done

for f in $OUTDIR/*.bg;
do
    bedGraphToBigWig "$f" dm6.chrom.sizes "$OUTDIR/$experiment.bw"
done;

#TF z-score
bamCoverage --bam ../processed_CR/ase-apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/ase-apkc_rep1.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 0.3204401
bamCoverage --bam ../processed_CR/ase-apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/ase-apkc_rep2.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 0.4509447
bamCoverage --bam ../processed_CR/btd_0hr_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/btd_0hr_rep1.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 7.158207
bamCoverage --bam ../processed_CR/btd_0hr_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/btd_0hr_rep2.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 6.506511
bamCoverage --bam ../processed_CR/FruC_0hr_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/FruC_0hr_rep1.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 0.1619853
bamCoverage --bam ../processed_CR/FruC_0hr_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/FruC_0hr_rep2.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 0.736927
bamCoverage --bam ../processed_CR/FruCOM_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/FruCOM_apkc_rep1.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 0.307137
bamCoverage --bam ../processed_CR/FruCOM_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/FruCOM_apkc_rep2.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 0.3530766
bamCoverage --bam ../processed_CR/Zld_0hr_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/Zld_0hr_rep1.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 0.7973307
bamCoverage --bam ../processed_CR/Zld_0hr_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/Zld_0hr_rep2.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 0.5361182
bamCoverage --bam ../processed_CR/Zld_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/Zld_apkc_rep1.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 1.177073
bamCoverage --bam ../processed_CR/Zld_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/Zld_apkc_rep2.small.bw --binSize 5 --maxFragmentLength 120 --numberOfProcessors 16 --scaleFactor 1.443169

#Histone Marks 

mkdir ../bigwigs_CR/TMM_norm
bamCoverage --bam ../processed_CR/H3K27ac_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/H3K27ac_apkc_rep1.bw --binSize 5 --scaleFactor 0.23342828 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K27ac_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/H3K27ac_apkc_rep2.bw --binSize 5 --scaleFactor 0.18624996 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K27me3_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/H3K27me3_apkc_rep1.bw --binSize 5 --scaleFactor 0.27023667 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K27me3_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/H3K27me3_apkc_rep2.bw --binSize 5 --scaleFactor 0.34675797 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K4me1_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/H3K4me1_apkc_rep1.bw --binSize 5 --scaleFactor 0.3348003 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K4me1_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/H3K4me1_apkc_rep2.bw --binSize 5 --scaleFactor 0.2346522 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K4me3_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/H3K4me3_apkc_rep1.bw --binSize 5 --scaleFactor 0.2003457 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K4me3_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/H3K4me3_apkc_rep2.bw --binSize 5 --scaleFactor 0.1739407 --numberOfProcessors 16

#Merge TMMs

INDIR="../bigwigs_CR/TMM_norm/"
OUTDIR="../bigwigs_CR/TMM_merged"

mkdir ../bigwigs_CR/TMM_merged
INDIR="../bigwigs_CR/TMM_norm"
OUTDIR="../bigwigs_CR/TMM_merged"
chmod 777 "$INDIR"
chmod 777 "$OUTDIR"
for rep1 in $INDIR/*_rep1.bw; do
  base=$(basename "$rep1")                 # e.g., H3K4me1_0hr_rep1.bw
  core="${base%_rep1.bw}"                  # e.g., H3K4me1_0hr
  rep2="$INDIR/${core}_rep2.bw"            # e.g., .../H3K4me1_0hr_rep2.bw
  out="$OUTDIR/${core}.bw"                 # e.g., bigwigs/TMM_merged/H3K4me1_0hr.bw

  if [[ -f "$rep2" ]]; then
	   echo "Merging: $base  +  $(basename "$rep2")  ->  $(basename "$out")"

		 bigwigCompare -b1 "$rep1" \
		 							 -b2 "$rep2" \
		 							 --operation mean --binSize 5 \
		 							 -p max -o "$out"
  else
    echo "WARNING: Missing pair for ${core}: ${rep2} not found. Skipping." >&2
  fi
done