# VLP-enrichment-methods-optimization-2025
Optimizing methods for virome analysis using a synthetic viral community
- ref_align.sh filters and maps the reads of all samples to all reference viral genomes from VirMock1, outputting .bam and .sam files. 
- vcontigs_align.sh filters and maps the reads of all the samples that have viral contigs (identified by Cenote-Taker2) to their respective viral contigs, outputting .bam and .sam files.
- coverage.py generates a table of mapped reads and total reads from .sam files

The following table provides the input files to run for each R script and the corresponding figure they generate. Input files can be found at https://zenodo.org/records/17436149.

| R script                                | Input file name                       | sheet name           | Figure generated   |
|:----------------------------------------|:--------------------------------------|:---------------------|:-------------------|
| Fig1_meta_VP.R                          | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result        | Fig1               |
| nan                                     | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage     | nan                |
| nan                                     | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result        | nan                |
| nan                                     | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage     | nan                |
| nan                                     | CHOPMC484-485_cenote_R_analysis.xlsx  | cenote_result        | nan                |
| nan                                     | CHOPMC484-485_cenote_R_analysis.xlsx  | viral_percentage     | nan                |
| Fig2_stool_spike-in_combined_analysis.R | metadata_250505.xlsx                  | nan                  | Fig2               |
| nan                                     | metadata_250617.xlsx                  | nan                  | nan                |
| nan                                     | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result        | nan                |
| nan                                     | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage     | nan                |
| nan                                     | 250505_nextseq_R_analysis.xlsx        | ref_viral_analysis   | nan                |
| nan                                     | 250505_nextseq_R_analysis.xlsx        | ref_viral_percentage | nan                |
| nan                                     | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result        | nan                |
| nan                                     | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage     | nan                |
| nan                                     | 250617_nextseq_R_analysis.xlsx        | ref_viral_analysis   | nan                |
| nan                                     | 250617_nextseq_R_analysis.xlsx        | ref_viral_percentage | nan                |
| Fig3_amplification.R                    | 250221_miniseq_ref_R_analysis.xlsx    | ref_viral_analysis   | Fig3               |
| nan                                     | 250221_miniseq_ref_R_analysis.xlsx    | ref_viral_percentage | nan                |
| nan                                     | 250228_nextseq_v2_ref_R_analysis.xlsx | ref_viral_analysis   | nan                |
| nan                                     | 250228_nextseq_v2_ref_R_analysis.xlsx | ref_viral_percentage | nan                |
| Fig4_nuclease.R                         | 250505_nextseq_R_analysis.xlsx        | ref_viral_analysis   | Fig4               |
| nan                                     | 250505_nextseq_R_analysis.xlsx        | ref_viral_percentage | nan                |
| nan                                     | 250516_nuc_titration_qpcr.xlsx        | nan                  | nan                |
| nan                                     | 250617_nextseq_R_analysis.xlsx        | ref_viral_analysis   | nan                |
| nan                                     | 250617_nextseq_R_analysis.xlsx        | ref_viral_percentage | nan                |
| Fig5A-B_300v1000.R                      | metadata_250505.xlsx                  | nan                  | Fig5A-B            |
| nan                                     | 300bp_vs_1000bp_R_analysis.xlsx       | ref_viral_analysis   | nan                |
| nan                                     | 300bp_vs_1000bp_R_analysis.xlsx       | ref_viral_percentage | nan                |
| nan                                     | 300bpvs1000bp_cenote_R_analysis.xlsx  | cenote_result        | nan                |
| nan                                     | 300bpvs1000bp_cenote_R_analysis.xlsx  | viral_percentage     | nan                |
| Fig5C-D_300vs1000len.R                  | metadata_250505.xlsx                  | nan                  | Fig5C-D            |
| nan                                     | 300bp_megahit_contig_table.tsv        | nan                  | nan                |
| nan                                     | 1000bp_megahit_contig_table.tsv       | nan                  | nan                |
| nan                                     | per_sample_bases.tsv                  | nan                  | nan                |
| nan                                     | 300bpvs1000bp_cenote_R_analysis.xlsx  | cenote_result        | nan                |
| Fig6_basemod.R                          | 250918_nextseq_R_analysis_v3.xlsx     | ref_viral_analysis   | Fig6               |
| nan                                     | 250918_nextseq_R_analysis_v3.xlsx     | ref_viral_percentage | nan                |
| nan                                     | 250922_T4+lambda_qpcr.xlsx            | nan                  | nan                |
| FigS1_meta_VP_Class.R                   | 250617_nextseq_cenote_R_analysis.xlsx | cenote_result        | Fig S1             |
| nan                                     | 250617_nextseq_cenote_R_analysis.xlsx | viral_percentage     | nan                |
| nan                                     | 250505_nextseq_cenote_R_analysis.xlsx | cenote_result        | nan                |
| nan                                     | 250505_nextseq_cenote_R_analysis.xlsx | viral_percentage     | nan                |
| nan                                     | CHOPMC484-485_cenote_R_analysis.xlsx  | cenote_result        | nan                |
| nan                                     | CHOPMC484-485_cenote_R_analysis.xlsx  | viral_percentage     | nan                |
| FigS2_BALsalivaOP.R                     | 250617_nextseq_R_analysis.xlsx        | ref_viral_analysis   | Fig S2             |
| nan                                     | 250617_nextseq_R_analysis.xlsx        | ref_viral_percentage | nan                |
| FigS3_LC-MS.R                           | LC-MS_data.xlsx                       | nan                  | Fig S3             |
