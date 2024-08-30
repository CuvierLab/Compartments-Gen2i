 #!/bin/bash
 
 # bedtools v2.29.1
 
 ## UNION PEAKSET
cat data/ChIPseq/hpl2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed data/ChIPseq/lin61_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed | bedtools sort -i - | bedtools merge   -i - >data/ChIPseq/hpl2_lin61_N2_narrow_union_peakset.bed 2> data/ChIPseq/hpl2_lin61_N2_narrow_union_peakset.log
cat data/ChIPseq/H3K9me2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed data/ChIPseq/H3K9me3_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed | bedtools sort -i - | bedtools merge   -i - >data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed 2> data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.log
cat data/ChIPseq/H3K9me2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed data/ChIPseq/H3K9me3_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed data/ChIPseq/lin61_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed data/ChIPseq/hpl2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed | bedtools sort -i - | bedtools merge   -i - >data/ChIPseq/H3K9me23_lin61_hpl2_optimal_N2_peakset.bed 2> data/ChIPseq/H3K9me23_lin61_hpl2_optimal_N2_peakset.log

 ## INTERSECT PEAKSET
 bedtools intersect   -a data/ChIPseq/lin61_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed -b data/ChIPseq/hpl2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed | cut -f1,2,3 >data/ChIPseq/hpl2_lin61_N2_sharp_intersect_peakset.bed 2> data/ChIPseq/hpl2_lin61_N2_sharp_intersect_peakset.log
