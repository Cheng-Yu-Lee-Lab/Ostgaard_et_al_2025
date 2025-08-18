#!/bin/sh
#SBATCH --job-name=job
#SBATCH --mail-user=rarjun@umich.edu
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4GB

module load python-anaconda3
source activate rarjun

#Make TF Peak Files

mkdir ../analysis_data/TF_peaks

cat "../../CUT_RUN_TOOLS/erm/5127-AR-15_GATAGCCA-CTGTGTTG_S54/peakcalling/macs2.narrow/5127-AR-15_GATAGCCA-CTGTGTTG_S54_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/erm/5127-AR-16_CCACAACA-TGTGGTAC_S55/peakcalling/macs2.narrow/5127-AR-16_CCACAACA-TGTGGTAC_S55_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/erm.bed

cat "../../CUT_RUN_TOOLS/ham/4899-AR-43_CTTCACTG-GTCCTAAG_S143/peakcalling/macs2.narrow/4899-AR-43_CTTCACTG-GTCCTAAG_S143_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/ham/4899-AR-44_CTCGACTT-GGTCAGAT_S144/peakcalling/macs2.narrow/4899-AR-44_CTCGACTT-GGTCAGAT_S144_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/ham.bed

cat "../../CUT_RUN_TOOLS/fruCOM_aPKC/6152-AR-10_GGAGGAAT-AACCGTTC_S33/peakcalling/macs2.narrow/6152-AR-10_GGAGGAAT-AACCGTTC_S33_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/fruCOM_aPKC/6152-AR-11_GACGTCAT-TCAACTGG_S34/peakcalling/macs2.narrow/6152-AR-11_GACGTCAT-TCAACTGG_S34_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/fruCOM_aPKC.bed

cat "../../CUT_RUN_TOOLS/zld_aPKC/6152-AR-3_TCATCTCC-CTAGCAAG_S26/peakcalling/macs2.narrow/6152-AR-3_TCATCTCC-CTAGCAAG_S26_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/zld_aPKC/6152-AR-4_CCAGTATC-ATCTCGCT_S27/peakcalling/macs2.narrow/6152-AR-4_CCAGTATC-ATCTCGCT_S27_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/zld_aPKC.bed

cat "../../CUT_RUN_TOOLS/fruCMYC/4362-AR-3_AGGAACAC-GACATTCC_S141/peakcalling/macs2.narrow/4362-AR-3_AGGAACAC-GACATTCC_S141_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/fruCMYC/5127-AR-27_GGTATAGG-CTGGAGTA_S65/peakcalling/macs2.narrow/5127-AR-27_GGTATAGG-CTGGAGTA_S65_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/fruC-myc.bed

cat "../../CUT_RUN_TOOLS/zld/EDL2_S27_L001/peakcalling/macs2.narrow/EDL2_S27_L001_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/zld/EDL4_S24_L001/peakcalling/macs2.narrow/EDL4_S24_L001_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/zld.bed

cat "../../CUT_RUN_TOOLS/pntP1/6044-AR-12_CAAGGTAC-GAGATACG_S107/peakcalling/macs2.narrow/6044-AR-12_CAAGGTAC-GAGATACG_S107_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/pntP1/6044-AR-13_AGACCTTG-GCACGTAA_S108/peakcalling/macs2.narrow/6044-AR-13_AGACCTTG-GCACGTAA_S108_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/pntP1.bed

cat "../../CUT_RUN_TOOLS/tll/5127-AR-11_GGAATGTC-GTGTTCCT_S52/peakcalling/macs2.narrow/5127-AR-11_GGAATGTC-GTGTTCCT_S52_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/tll/5127-AR-12_TGGTGAAG-GCTGTAAG_S53/peakcalling/macs2.narrow/5127-AR-12_TGGTGAAG-GCTGTAAG_S53_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/tll.bed

cat "../../CUT_RUN_TOOLS/btd/6044-AR-14_GTCGTTAC-GCTTAGCT_S109/peakcalling/macs2.narrow/6044-AR-14_GTCGTTAC-GCTTAGCT_S109_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/btd/6044-AR-15_GTAACCGA-GGTGTCTT_S110/peakcalling/macs2.narrow/6044-AR-15_GTAACCGA-GGTGTCTT_S110_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/btd.bed

cat "../../CUT_RUN_TOOLS/ase_aPKC/10379-CO-13_S13/peakcalling/macs2.narrow/10379-CO-13_S13_peaks.narrowPeak" \
"../../CUT_RUN_TOOLS/ase_aPKC/10379-CO-13_S13/peakcalling/macs2.narrow/10379-CO-13_S13_peaks.narrowPeak" \
| bedtools sort -i stdin | bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v | bedtools merge -i stdin > ../analysis_data/TF_peaks/ase.bed

mkdir ../analysis_data/flowchart_v4

grep -v 'gene_id .*RNA:' dm6.ncbiRefSeq.gtf > filtered_dm6.ncbiRefSeq.gtf

awk 'BEGIN{OFS="\t"} $3=="transcript" {
    if ($7=="+") {
        start = ($4 - 250 > 0) ? $4 - 250 : 1; 
        end = ($4 + 250);
    } else {
        start = ($5 - 250 > 0) ? $5 - 250 : 1;
        end = ($5 + 250);
    }
    print $1, start, end, $9, ".", $7;
}' filtered_dm6.ncbiRefSeq.gtf > ../analysis_data/flowchart_v4/dm6.promoters.bed

#Top level of Flowchart
bedtools intersect -a "../peaks_ATAC/0hr_rep1/0hr_rep1_peaks.narrowPeak" \
	-b "../peaks_ATAC/0hr_rep2/0hr_rep2_peaks.narrowPeak" -wa -u \
	| bedtools intersect -a stdin -b "../peaks_ATAC/0hr_rep3/0hr_rep3_peaks.narrowPeak" -wa -u \
	| bedtools intersect -a stdin -b dm6.blacklist.noCHR.bed -v \
	| bedtools intersect -a stdin -b ../analysis_data/flowchart_v4/dm6.promoters.bed -v  \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/typeII-promoter.bed

cat ../analysis_data/flowchart_v4/typeII-promoter.bed | \
	   awk -F'/t' -v OFS='/t' '{ $1 = "chr" $1 }1' | \
	   annotatePeaks.pl - dm6 > ../analysis_data/flowchart_v4/typeII-promoter.anno.txt

#Next level of flowchart (TF Intersection)
bedtools intersect -a "../analysis_data/flowchart_v4/typeII-promoter.bed" \
	-b "../analysis_data/TF_peaks/fruC-myc.bed" -wa -u \
	| bedtools intersect -a stdin -b "../analysis_data/TF_peaks/zld.bed" -wa -u \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/regulatory-tf_t2.bed

#Next level of flowchart (INP Intersection)
bedtools intersect -a "../analysis_data/flowchart_v4/typeII-promoter.bed" \
	-b "../peaks_ATAC/24hr_rep1/24hr_rep1_peaks.narrowPeak" -wa -u \
	| bedtools intersect -a stdin -b "../peaks_ATAC/24hr_rep2/24hr_rep2_peaks.narrowPeak" -wa -u \
	| bedtools intersect -a stdin -b "../peaks_ATAC/24hr_rep3/24hr_rep3_peaks.narrowPeak" -wa -u \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/INP_open.bed

bedtools intersect -a "../analysis_data/flowchart_v4/regulatory-tf_t2.bed" \
	-b "../analysis_data/flowchart_v4/INP_open.bed" -wa -u \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/regulatory-tf_t2_INP_open.bed

bedtools intersect -a "../analysis_data/flowchart_v4/regulatory-tf_t2.bed" \
	-b "../analysis_data/flowchart_v4/INP_open.bed" -wa -v \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/regulatory-tf_t2_INP_closed.bed

#Type I 
bedtools intersect -a "../analysis_data/flowchart_v4/typeII-promoter.bed" \
	-b "../peaks_ATAC/apkc_rep1/apkc_rep1_peaks.narrowPeak" -wa -u \
	| bedtools intersect -a stdin -b "../peaks_ATAC/apkc_rep2/apkc_rep2_peaks.narrowPeak" -wa -u \
	| bedtools intersect -a stdin -b "../peaks_ATAC/apkc_rep3/apkc_rep3_peaks.narrowPeak" -wa -u \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/typeI_open.bed

#Type I (right side)
bedtools intersect -a "../analysis_data/flowchart_v4/regulatory-tf_t2_INP_closed.bed" \
	-b "../analysis_data/flowchart_v4/typeI_open.bed" -wa -u \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/Pan-NB.bed

bedtools intersect -a "../analysis_data/flowchart_v4/regulatory-tf_t2_INP_closed.bed" \
	-b "../analysis_data/flowchart_v4/typeI_open.bed" -wa -v \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/T2_NB.bed

#Type I (left side)
bedtools intersect -a "../analysis_data/flowchart_v4/typeI_open.bed" \
	-b "../analysis_data/TF_peaks/fruCOM_aPKC.bed" -wa -u \
	| bedtools intersect -a stdin -b "../analysis_data/TF_peaks/zld_aPKC.bed" -wa -u \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/regulatory-tf_t1.bed

bedtools intersect -a "../analysis_data/flowchart_v4/regulatory-tf_t2_INP_open.bed" \
	-b "../analysis_data/flowchart_v4/regulatory-tf_t1.bed" -wa -u \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/Pan-neural.bed

bedtools intersect -a "../analysis_data/flowchart_v4/regulatory-tf_t2_INP_open.bed" \
	-b "../analysis_data/flowchart_v4/regulatory-tf_t1.bed" -wa -v \
	| bedtools sort -i stdin \
	| bedtools merge -i stdin > ../analysis_data/flowchart_v4/T2_lineage.bed