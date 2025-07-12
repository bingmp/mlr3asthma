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

## Mlr3

The folder named mlr3 is the main bit of the dataset federated data
analysis process, including machine learning such as random forests and
logistics regression, as well as shap analysis and more.

    list.files('3_mlr3/')

    ## [1] "code"    "rawdata" "RDS"     "result"

    list.files('3_mlr3/code/')

    ## [1] "step_01_upset.R"       "step_02_coDEGS_exp.R"  "step_03_5fold_cross.R"
    ## [4] "step_04_train_test.R"  "step_05_machlearn.R"   "step_06_logis.R"      
    ## [7] "step_07_DCA.R"         "step_08_rf_shap.R"

## Packages infomation

    sessionInfo()

    ## R version 4.5.1 (2025-06-13)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=zh_CN.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=zh_CN.UTF-8        LC_COLLATE=zh_CN.UTF-8    
    ##  [5] LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=zh_CN.UTF-8   
    ##  [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Asia/Shanghai
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.5.1    fastmap_1.2.0     cli_3.6.5         tools_4.5.1      
    ##  [5] htmltools_0.5.8.1 rstudioapi_0.17.1 yaml_2.3.10       rmarkdown_2.29   
    ##  [9] knitr_1.50        xfun_0.52         digest_0.6.37     rlang_1.1.6      
    ## [13] evaluate_1.0.4
