# GTEX_liver_WGCNA_Human


three files in Human GTEX liver WGCNA
(1)  run locally to get working image
(2)  run in cluster to get modules
(3)  run locally to get GO enrichment for each module colors  



1. run Gtex_local with gtex.txt file to get working image saved with all data (network_expr_data.RData)
2. run in cluster (using wg_gtex.sub and wg_gtex.sh)  using Gtex_combined_R together with previously saved working RData
3. download the 0-gtex, 1-gtex and 2-gtex.RData files together with the consensusAnalysis.csv file (with module information)
4. run locally with GTEX_Go with previously downloaded working RData file and using Annotation.csv file as reference.
5. get the GO.csv file with all GO module with ModuleColors.
