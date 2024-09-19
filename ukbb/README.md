# Code to reproduce UK Biobank analysis

Overview of files:

- `01-preprocessing.R`: preprocesses UK Biobank data 
- `02-subset.R`: downsamples UK Biobank data
- `03-merge.R`: combines UK Biobank significance results
- `04-LD-blocks.R`: estimate independent LD blocks
- `05-sffdr.R`: applying sfFDR to GWAS of BMI with obesity-related informative traits
- `06-sffdr-null.R`: applying sfFDR to GWAS of BMI with non-informative traits (permuted) in the UK Biobank data 
    
The above code uses `../00-helper.R`.
