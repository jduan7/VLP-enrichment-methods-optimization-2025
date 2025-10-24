# VLP-enrichment-methods-optimization-2025
Optimizing methods for virome analysis using a synthetic viral community
- ref_align.sh filters and maps the reads of all samples to all reference viral genomes from VirMock1, outputting .bam and .sam files. 
- vcontigs_align.sh filters and maps the reads of all the samples that have viral contigs (identified by Cenote-Taker2) to their respective viral contigs, outputting .bam and .sam files.
- coverage.py generates a table of mapped reads and total reads from .sam files

The following table provides the input files to run for each R script and the corresponding figure they generate. Input files can be found at https://zenodo.org/records/17436149.

| R script                                | Input file name                       | sheet name           | Figure generated   |
|:----------------------------------------|:--------------------------------------|:---------------------|:-------------------|
| Fig1_meta_VP.R                          | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result        | Fig1               |
|                                         | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage     |                    |
|                                         | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result        |                    |
|                                         | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage     |                    |
|                                         | CHOPMC484-485_cenote_R_analysis.xlsx  | cenote_result        |                    |
|                                         | CHOPMC484-485_cenote_R_analysis.xlsx  | viral_percentage     |                    |
| Fig2_stool_spike-in_combined_analysis.R | metadata_250505.xlsx                  | NA                   | Fig2               |
|                                         | metadata_250617.xlsx                  | NA                   |                    |
|                                         | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result        |                    |
|                                         | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage     |                    |
|                                         | 250505_nextseq_R_analysis.xlsx        | ref_viral_analysis   |                    |
|                                         | 250505_nextseq_R_analysis.xlsx        | ref_viral_percentage |                    |
|                                         | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result        |                    |
|                                         | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage     |                    |
|                                         | 250617_nextseq_R_analysis.xlsx        | ref_viral_analysis   |                    |
|                                         | 250617_nextseq_R_analysis.xlsx        | ref_viral_percentage |                    |
| Fig3_amplification.R                    | 250221_miniseq_ref_R_analysis.xlsx    | ref_viral_analysis   | Fig3               |
|                                         | 250221_miniseq_ref_R_analysis.xlsx    | ref_viral_percentage |                    |
|                                         | 250228_nextseq_v2_ref_R_analysis.xlsx | ref_viral_analysis   |                    |
|                                         | 250228_nextseq_v2_ref_R_analysis.xlsx | ref_viral_percentage |                    |
| Fig4_nuclease.R                         | 250505_nextseq_R_analysis.xlsx        | ref_viral_analysis   | Fig4               |
|                                         | 250505_nextseq_R_analysis.xlsx        | ref_viral_percentage |                    |
|                                         | 250516_nuc_titration_qpcr.xlsx        | NA                   |                    |
|                                         | 250617_nextseq_R_analysis.xlsx        | ref_viral_analysis   |                    |
|                                         | 250617_nextseq_R_analysis.xlsx        | ref_viral_percentage |                    |
| Fig5A-B_300v1000.R                      | metadata_250505.xlsx                  | NA                   | Fig5A-B            |
|                                         | 300bp_vs_1000bp_R_analysis.xlsx       | ref_viral_analysis   |                    |
|                                         | 300bp_vs_1000bp_R_analysis.xlsx       | ref_viral_percentage |                    |
|                                         | 300bpvs1000bp_cenote_R_analysis.xlsx  | cenote_result        |                    |
|                                         | 300bpvs1000bp_cenote_R_analysis.xlsx  | viral_percentage     |                    |
| Fig5C-D_300vs1000len.R                  | metadata_250505.xlsx                  | NA                   | Fig5C-D            |
|                                         | 300bp_megahit_contig_table.tsv        | NA                   |                    |
|                                         | 1000bp_megahit_contig_table.tsv       | NA                   |                    |
|                                         | per_sample_bases.tsv                  | NA                   |                    |
|                                         | 300bpvs1000bp_cenote_R_analysis.xlsx  | cenote_result        |                    |
| Fig6_basemod.R                          | 250918_nextseq_R_analysis_v3.xlsx     | ref_viral_analysis   | Fig6               |
|                                         | 250918_nextseq_R_analysis_v3.xlsx     | ref_viral_percentage |                    |
|                                         | 250922_T4+lambda_qpcr.xlsx            | NA                   |                    |
| FigS1_meta_VP_Class.R                   | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result        | Fig S1             |
|                                         | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage     |                    |
|                                         | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result        |                    |
|                                         | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage     |                    |
|                                         | CHOPMC484-485_cenote_R_analysis.xlsx  | cenote_result        |                    |
|                                         | CHOPMC484-485_cenote_R_analysis.xlsx  | viral_percentage     |                    |
| FigS2_BALsalivaOP.R                     | 250617_nextseq_R_analysis.xlsx        | ref_viral_analysis   | Fig S2             |
|                                         | 250617_nextseq_R_analysis.xlsx        | ref_viral_percentage |                    |
| FigS3_LC-MS.R                           | LC-MS_data.xlsx                       | NA                   | Fig S3             |
