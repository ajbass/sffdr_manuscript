# Code to reproduce UK Biobank analysis

Overview of files:

    - `01-preprocessing.R`: preprocesses UK Biobank data 
    - `03-subset.R`: downsamples UK Biobank data
    - `04-merge.R`: combines UK Biobank significance results
    - `05-LD-blocks.R`: estimate independent LD blocks
    - `06-sffdr.R`: applying sfFDR to GWAS of BMI with obesity-related informative traits
    - `07-sffdr-null.R`: applying sfFDR to GWAS of BMI with non-informative traits (permuted) in the UK Biobank data 
    
The above code uses `../00-helper.R`. 
