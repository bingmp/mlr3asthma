## GSE152004

The folder named GSE152004 is the main bit of the data analysis process
for the dataset GSE152004, including WGCNA analysis, differential gene
analysis, volcano and heat maps, and more.

    list.files('1_GSE152004/')

    ## [1] "code"    "rawdata" "RDS"     "result"

    list.files('1_GSE152004/code/')

    ##  [1] "step_00_packages.R"                "step_01_data.R"                   
    ##  [3] "step_02_vstMat_rm.R"               "step_03_pickSoftThreshold.R"      
    ##  [5] "step_04_WGCNA.R"                   "step_05_hclust.R"                 
    ##  [7] "step_06_DEseq2.R"                  "step_07_DEseq2_DEGs.R"            
    ##  [9] "step_08_DEseq2_plot_boxplot.R"     "step_08_DEseq2_plot_correlation.R"
    ## [11] "step_08_DEseq2_plot_heatmap.R"     "step_08_DEseq2_plot_PCA.R"        
    ## [13] "step_08_DEseq2_plot_ROC.R"         "step_08_DEseq2_plot_volcano.R"

## GSE67472

The folder named GSE67472 is the main bit of the data analysis process
for the dataset GSE67472, including differential gene analysis, volcano
and heat maps, and more.

    list.files('2_GSE67472/')

    ## [1] "code"    "rawdata" "RDS"     "result"

    list.files('2_GSE67472/code/')

    ## [1] "step_01_data.R"            "step_02_geneID.R"         
    ## [3] "step_03_norm.R"            "step_04_DEGs.R"           
    ## [5] "step_05_plot_correlatin.R" "step_05_plot_heatmap.R"   
    ## [7] "step_05_plot_PCA.R"        "step_05_plot_ROC.R"       
    ## [9] "step_05_plot_volcano.R"

## mlr3

The folder named mlr3 is the main bit of the dataset federated data
analysis process, including machine learning such as random forests and
logistics regression, as well as shap analysis and more.

    list.files('3_mlr3/')

    ## [1] "code"    "rawdata" "RDS"     "result"

    list.files('3_mlr3/code/')

    ## [1] "step_01_upset.R"       "step_02_coDEGS_exp.R"  "step_03_5fold_cross.R"
    ## [4] "step_04_train_test.R"  "step_05_machlearn.R"   "step_06_logis.R"      
    ## [7] "step_07_DCA.R"         "step_08_rf_shap.R"
