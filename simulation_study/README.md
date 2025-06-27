# Code to reproduce simulation study

Overview of files:

- `00-generate_data.R`: functions to generate data
- `01-sffdr-independent.R`: run independent SNP setting for sfFDR
- `02-fdr-methods-independent.R`: run independent SNP setting for AdaPT, CAMT, and Boca-Leek
- `03-sffdr-independent-null.R`: run null independent SNP setting for sfFDR
- `04-fdr-methods-independent-null.R`: run null independent SNP setting for AdaPT, CAMT, and Boca-Leek
- `05-sffdr-dependent.R`: run dependent SNP setting for sfFDR
- `06-computational-time.R`: run computational time comparisons for sfFDR

The above code uses `../00-helper.R`. 
