# TCGA_GDSC_drug_response_modeling

Pre-processing of TCGA RNA-seq expression data or GDSC Microarray expression data to build classification models (ML/DL) for drug response in different cancers.

For RNA-seq data from TCGA:
Step 1: Run process_RNAseq_DEgenes_ML_100cuts.R - this generates 100 cuts of scaled train and test data sets based on sig. DE genes, and building multiple machine learning models.
Step 2: Run DNN_h2o_binary_classf.R - this builds DL models using H2O framework on the scaled data sets generated from step 1.

For Microarray data from GDSC:
Step 1: Run MA.end.to.end.proc.R - this generates 100 cuts of scaled train and test data sets based on sig. DE genes.
Step 2a: Run run_ML_multirun_affy.R - this builds multiple ML models on each cut of train and test generated from step 1.
Step 2b: Run DNN_h2o_binary_classf.R - this builds DL models using H2O framework on the scaled data sets generated from step 1.

Note: link to raw and processed data sets:
https://appliedbiomath.box.com/s/0lgnc93iq47xetrwspcuwnx64f4qe1pk
please contact sgujja@appliedbiomath.com for any questions.
