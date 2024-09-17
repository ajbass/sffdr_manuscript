# Code to reproduce EGPA analysis

Overview of files:

    - `01-preprocessing.R`: preprocesses EGPA + informative studies
    - `02-sffdr-egpa.R`: apply sffdr to EGPA w/ EGPA-informative traits
    - `03-sffdr-null-permuted`: apply sffdr to EGPA w/ null (permuted) obesity-related traits
    - `04-sffdr-null-original`: apply sffdr to EGPA w/ null (unpermuted) obesity-related traits

The above code uses `../00-helper.R`. 
