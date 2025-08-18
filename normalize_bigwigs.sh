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

mkdir ../bigwigs_ATAC/z-score
mkdir ../bigwigs_ATAC/z-score/bg

for f in $OUTDIR/*.bg;
do
    experiment=$(basename $f | awk -F '.' '{print $1}')

    mean=$(awk '{ total += $4 } END { print total/NR }' "$f")
    std=$(awk -v mean="$mean" '{ total += ($4 - mean)^2 } END { print sqrt(total/NR) }' "$f")
    awk -v mean="$mean" -v std="$std" '{ printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, ($4 - mean) / std }' "$f" > "../bigwigs_ATAC/z-score/bg/$experiment.z-score.bg"
    bedGraphToBigWig "../bigwigs_ATAC/z-score/bg/$experiment.z-score.bg" dm6.chrom.sizes "../bigwigs_ATAC/z-score/$experiment.z-score.bw"
done;

#TF z-score

mkdir ../bigwigs_CR/TF_merged
INDIR="../bigwigs_CR/"
OUTDIR="../bigwigs_CR/TF_merged"

for rep1 in $INDIR/*_rep1.small.bw; do
  base=$(basename "$rep1")                
  core="${base%_rep1.small.bw}"                
  rep2="$INDIR/${core}_rep2.small.bw"           
  out="$OUTDIR/${core}.bg"                

  if [[ -f "$rep2" ]]; then
    echo "Merging: $base  +  $(basename "$rep2")  ->  $(basename "$out")"
	wiggletools write_bg "$out" mean \
                 "$rep1" \
                 "$rep2" 
  else
    echo "WARNING: Missing pair for ${core}: ${rep2} not found. Skipping." >&2
  fi
done

mkdir ../bigwigs_CR/z-score
mkdir ../bigwigs_CR/z-score/bg

for f in $OUTDIR/*.bg;
do
    experiment=$(basename $f | awk -F '.' '{print $1}')

    mean=$(awk '{ total += $4 } END { print total/NR }' "$f")
    std=$(awk -v mean="$mean" '{ total += ($4 - mean)^2 } END { print sqrt(total/NR) }' "$f")
    awk -v mean="$mean" -v std="$std" '{ printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, ($4 - mean) / std }' "$f" > "../bigwigs_CR/z-score/bg/$experiment.z-score.bg"
    bedGraphToBigWig "../bigwigs_CR/z-score/bg/$experiment.z-score.bg" dm6.chrom.sizes "../bigwigs_CR/z-score/$experiment.z-score.bw"
done;

mkdir ../bigwigs_CR/TMM_norm
bamCoverage --bam ../processed_CR/H3K27ac_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/H3K27ac_apkc_rep1.bw --binSize 5 --scaleFactor 0.23342828 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K27ac_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/H3K27ac_apkc_rep2.bw --binSize 5 --scaleFactor 0.18624996 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K27me3_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/H3K27me3_apkc_rep1.bw --binSize 5 --scaleFactor 0.27023667 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K27me3_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/H3K27me3_apkc_rep2.bw --binSize 5 --scaleFactor 0.34675797 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K4me1_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/H3K4me1_apkc_rep1.bw --binSize 5 --scaleFactor 0.3348003 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K4me1_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/H3K4me1_apkc_rep2.bw --binSize 5 --scaleFactor 0.2346522 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K4me3_apkc_rep1/*.bam  -o ../bigwigs_CR/TMM_norm/H3K4me3_apkc_rep1.bw --binSize 5 --scaleFactor 0.2003457 --numberOfProcessors 16
bamCoverage --bam ../processed_CR/H3K4me3_apkc_rep2/*.bam  -o ../bigwigs_CR/TMM_norm/H3K4me3_apkc_rep2.bw --binSize 5 --scaleFactor 0.1739407 --numberOfProcessors 16

mkdir ../bigwigs_CR/TMM_merged
INDIR="../bigwigs_CR/TMM_norm/"
OUTDIR="../bigwigs_CR/TMM_merged"

for rep1 in $INDIR/H3*_rep1.bw; do
  base=$(basename "$rep1")              
  core="${base%_rep1.bw}"               
  rep2="$INDIR/${core}_rep2.bw"           
  out="$OUTDIR/${core}.bg"                

  if [[ -f "$rep2" ]]; then
    echo "Merging: $base  +  $(basename "$rep2")  ->  $(basename "$out")"
	wiggletools write_bg "$out" mean \
                 "$rep1" \
                 "$rep2" 
  else
    echo "WARNING: Missing pair for ${core}: ${rep2} not found. Skipping." >&2
  fi
done

mkdir ../bigwigs_CR/z-score
mkdir ../bigwigs_CR/z-score/bg

for f in $OUTDIR/H3*.bg;
do
    experiment=$(basename $f | awk -F '.' '{print $1}')

    mean=$(awk '{ total += $4 } END { print total/NR }' "$f")
    std=$(awk -v mean="$mean" '{ total += ($4 - mean)^2 } END { print sqrt(total/NR) }' "$f")
    awk -v mean="$mean" -v std="$std" '{ printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, ($4 - mean) / std }' "$f" > "../bigwigs_CR/z-score/bg/$experiment.z-score.bg"
    bedGraphToBigWig "../bigwigs_CR/z-score/bg/$experiment.z-score.bg" dm6.chrom.sizes "../bigwigs_CR/z-score/$experiment.z-score.bw"
done;