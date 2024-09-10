 #!/bin/bash
 
 # deeptools 3.5.1
 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/H3K9me2_hpl2.lin61_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_hpl2.lin61_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_hpl2.lin61_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/hpl2_hpl2.lin61_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/hpl2_hpl2.lin61_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/hpl2_hpl2.lin61_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/LEM2_hpl2.lin61_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/LEM2_hpl2.lin61_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/LEM2_hpl2.lin61_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/H3K9me2_I158A_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_I158A_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_I158A_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/lin61_I158A_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/lin61_I158A_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/lin61_I158A_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/hpl2_I158A_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/hpl2_I158A_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/hpl2_I158A_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/LEM2_I158A_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/LEM2_I158A_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/LEM2_I158A_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/H3K9me2_met2.set25.set32_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_met2.set25.set32_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_met2.set25.set32_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/lin61_met2.set25.set32_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/lin61_met2.set25.set32_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/lin61_met2.set25.set32_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/hpl2_met2.set25.set32_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/hpl2_met2.set25.set32_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/hpl2_met2.set25.set32_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/H3K9me2_hpl2_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_hpl2_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_hpl2_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/lin61_hpl2_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/lin61_hpl2_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/lin61_hpl2_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/LEM2_hpl2_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/LEM2_hpl2_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/LEM2_hpl2_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/H3K9me2_lin61_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_lin61_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/H3K9me2_lin61_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/lin61_lin61_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/lin61_lin61_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/lin61_lin61_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/hpl2_lin61_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/hpl2_lin61_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/hpl2_lin61_newzs_prof.matrix.tab 
computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/LEM2_lin61_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/LEM2_lin61_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/LEM2_lin61_newzs_prof.matrix.tab
