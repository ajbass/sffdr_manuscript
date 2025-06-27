# Code to reproduce UK Biobank analysis

Overview of files:

- `01-preprocessing.R`: preprocesses UK Biobank data 
- `02-subset.R`: downsamples UK Biobank data
- `03-filtering.R`: combines UK Biobank significance results
- `04-LD-blocks.R`: estimate independent LD blocks
- `05-bmi-inform-sampling.R`: applying sfFDR to GWAS of BMI with obesity-related informative traits (informative sampling)
- `06-bmi-null-case1.R`: applying sfFDR to GWAS of BMI with non-informative traits (permuted)
- `06b-bmi-null-case2.R`: applying sfFDR to GWAS of permuted BMI with informative traits
- `06c-bmi-null-case3.R`: applying sfFDR to GWAS of BMI (a mixture of null/non-null SNPs) with informative traits  
- `07-bmi-mixture-null-non-null.R`: applying sfFDR to GWAS of BMI with a mixture of null/non-null informative traits 
- `08-finngen.R`: applying sfFDR to GWAS of BMI in the UK Biobank using BMI from FinnGen as an informative trait
- `09-giant.R`: applying sfFDR to GWAS of BMI in the UK Biobank using BMI from GIANT as an informative trait
- `10-mvp.R`: applying sfFDR to GWAS of BMI in the UK Biobank using BMI from MVP as an informative trait
- `11-prune.R`/`11b-prune-null.R`/`11c-prune-null.R`: Evaluating sfFDR using pruned SNPs to train the model
    
The above code uses `../00-helper.R`.
