# Code to reproduce EGPA analysis

Overview of files:

- `01-preprocessing.R`: preprocesses EGPA + informative studies
- `02-egpa.R`: apply sfFDR to EGPA w/ EGPA-informative traits
- `02b-egpa-null-permuted-traits.R`: apply sfFDR to EGPA w/ null (permuted) obesity-related traits
- `02c-egpa-null-original-traits.R`: apply sfFDR to EGPA w/ null (unpermuted) obesity-related traits
- `03-myositis-finngen.R`: apply sfFDR to GWAS of myositis (FinnGen)
- `04-juven-finngen.R`: apply sfFDR to GWAS of juvenile arthritis (FinnGen)
- `05-sle-finngen.R`: apply sfFDR to GWAS of systemic lupus erythematosus (FinnGen)
- `06-autoimmune-thyroiditis-finngen.R`: apply sfFDR to GWAS of autoimmune thyroiditis (FinnGen)

The above code uses `../00-helper.R`. 
