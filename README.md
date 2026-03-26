# VLP-enrichment-methods-optimization-2025
Optimizing methods for virome analysis using a synthetic viral community
- ref_align.sh filters and maps the reads of all samples to all reference viral genomes from VirMock1, outputting .bam and .sam files. 
- vcontigs_align.sh filters and maps the reads of all the samples that have viral contigs (identified by Cenote-Taker2) to their respective viral contigs, outputting .bam and .sam files.
- coverage.py generates a table of mapped reads and total reads from .sam files

The following table provides the input files to run for each R script and the corresponding figure they generate. Input files can be found at https://zenodo.org/records/18750187.

## Figure Scripts and Input Files

| R script | Input file name | Sheet name | Figure |
|----------|-----------------|------------|--------|
| Fig1_FigS1_meta_VP.R | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result | Fig 1, Fig S1 |
|  | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | CHOPMC484-485_cenote_R_analysis.xlsx | cenote_result |  |
|  | CHOPMC484-485_cenote_R_analysis.xlsx | viral_percentage |  |
| Fig2_FigS5_stoolVP1-3.R | metadata_250505 |  | Fig 2, Fig S5 |
|  | metadata_250617 |  |  |
|  | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | 250505_nextseq_R_analysis.xlsx | ref_viral_analysis |  |
|  | 250505_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
|  | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_analysis |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
| Fig3_amplification.R | 250221_miniseq_ref_R_analysis.xlsx | ref_viral_analysis | Fig 3 |
|  | 250221_miniseq_ref_R_analysis.xlsx | ref_viral_percentage |  |
|  | 250228_nextseq_v2_ref_R_analysis.xlsx | ref_viral_analysis |  |
|  | 250228_nextseq_v2_ref_R_analysis.xlsx | ref_viral_percentage |  |
| Fig4_nuclease.R | 250505_nextseq_R_analysis.xlsx | ref_viral_analysis | Fig 4 |
|  | 250505_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
|  | 250516_nuc_titration_qpcr.xlsx |  |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_analysis |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
| Fig5A-F_300v1000.R | metadata_250505 |  | Fig 5A-F |
|  | 300bp_vs_1000bp_R_analysis.xlsx | ref_viral_analysis |  |
|  | 300bp_vs_1000bp_R_analysis.xlsx | ref_viral_percentage |  |
|  | 300bpvs1000bp_cenote_R_analysis.xlsx | cenote_result |  |
|  | 300bpvs1000bp_cenote_R_analysis.xlsx | viral_percentage |  |
| Fig5G-H_300vs1000len.R | metadata_250505 |  | Fig 5G-H |
|  | 300bp_megahit_contig_table.tsv |  |  |
|  | 1000bp_megahit_contig_table.tsv |  |  |
|  | per_sample_bases.tsv |  |  |
|  | 300bpvs1000bp_cenote_R_analysis.xlsx | cenote_result |  |
| Fig6_basemod.R | 250918_nextseq_R_analysis_v3.xlsx | ref_viral_analysis | Fig 6 |
|  | 250918_nextseq_R_analysis_v3.xlsx | ref_viral_percentage |  |
|  | 250922_T4+lambda_qpcr.xlsx |  |  |
| FigS2_meta_VP_geNomad.R | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result | Fig S2 |
|  | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | CHOPMC484-485_cenote_R_analysis.xlsx | cenote_result |  |
|  | CHOPMC484-485_cenote_R_analysis.xlsx | viral_percentage |  |
| FigS3_spiked-in_BALsalivaOP.R | metadata_250617 |  | Fig S3 |
|  | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_analysis |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
| FigS4_stool_spike-in_pcoa.R | metadata_250505 |  | Fig S4 |
|  | metadata_250617 |  |  |
|  | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | 250505_nextseq_R_analysis.xlsx | ref_viral_analysis |  |
|  | 250505_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
|  | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_analysis |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
| FigS6_KrakenResult.Rmd | samples_for_kraken.xlsx |  | Fig S6 |
|  | all_samples.tsv |  |  |
| FigS7_BALsalivaOP.R | metadata_250617 |  | Fig S7 |
|  | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_analysis |  |
|  | 250617_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
| FigS8_2ndstrand.R | 260116_nextseq_R_analysis.xlsx | ref_viral_analysis | Fig S8A-N |
|  | 260116_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
|  | 260116_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 260116_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
| FigS9_2ndstrand_usedist_stats.R | 260116_nextseq_R_analysis.xlsx | ref_viral_analysis | Fig S8O |
|  | 260116_nextseq_R_analysis.xlsx | ref_viral_percentage |  |
|  | 260116_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 260116_nextseq_cenote_R_analysis.xlsx | viral_percentage |  |
| FigS10_seq_depth_viral_contig.R | 250221_miniseq_cenote_R_analysis.xlsx | cenote_result | Fig S9 |
|  | 250228_nextseq_v2_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250815_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 250918_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | 260116_nextseq_cenote_R_analysis.xlsx | cenote_result |  |
|  | Revised_SuppTable_v8.xlsx | S15_sample_metadata |  |
| FigS9_LC-MS.R | LC-MS_data.xlsx |  | Fig S9 |
