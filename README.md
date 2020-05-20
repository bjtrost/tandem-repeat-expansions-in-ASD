# tandem-repeat-expansions-in-ASD

The following is a list of scripts meant to be invoked directly (as opposed to being called by other scripts):

- generate_EH_genotype_table.py: Generate a table containing ExpansionHunter genotypes given multiple ExpansionHunter VCF files as input.
- EHDN_raw_call_count.py: Output the raw number of calls (anchored in-repeat reads) from an ExpansionHunter Denovo output file.
- EHDN_call_count_boxplots.R: Generate boxplots summarizing the call counts, stratified by DNA library preparation method and sequencing platform and by ethnicity.
- EHDN_STR_summary.sh: Generate plots summarizing the distribution of repeat units (motifs) detected by ExpansionHunter Denovo.

The following is a list of scripts that are unlikely to be useful to the user, but that document the exact procedure used to perform an analysis in the paper.

- find_EH_EHDN_correlation.py: Documents the procedure used to determine the correlation between ExpansionHunter and ExpansionHunter Denovo calls
- EHDN_validation_TRF_overlap.py: Documents the procedure used to determine the validation rate of ExpansionHunter Denovo calls based on overlap with TRF calls and Asmvar calls (derived from PacBio data) in the HuRef genome.

The following is a list of library files, which contain helper functions used by other scripts in this repository:

- BTlib.py (Python functions)
- BTlib/R/BTlib.R (R functions)

