# tandem-repeat-expansions-in-ASD
This repository contains code used to perform analyses done in the paper "Genome-wide detection of tandem DNA repeats expanded in autism".

The following is a list of scripts meant to be invoked directly (as opposed to being called by other scripts):

- generate_EH_genotype_table.py: Generate a table containing ExpansionHunter genotypes given multiple ExpansionHunter VCF files as input.
- EHDN_raw_call_count.py: Output the raw number of calls (anchored in-repeat reads) from an ExpansionHunter Denovo output file.
- EHDN_call_count_boxplots.R: Generate boxplots summarizing the call counts, stratified by DNA library preparation method and sequencing platform and by ethnicity.
- EHDN_STR_summary.sh: Generate plots summarizing the distribution of repeat units (motifs) detected by ExpansionHunter Denovo.
- DBSCAN.EHdn.parallel.R: Identify expansions from EHdn genotype table


The following is a list of scripts that are unlikely to be useful to the user, but that document the exact procedure used to perform an analysis in the paper.

- find_EH_EHDN_correlation.py: Documents the procedure used to determine the correlation between ExpansionHunter and ExpansionHunter Denovo calls
- EHDN_validation_TRF_overlap.py: Documents the procedure used to determine the validation rate of ExpansionHunter Denovo calls based on overlap with TRF calls and Asmvar calls (derived from PacBio data) in the HuRef genome.
- denovoSNVCNVAnalysis.R: Test the association between expansions and other rare loss-of-function variants.
- disease.loci.test.R: Fisher's Exact test for disease loci burden comparison between affected and unaffected individuals
- genPlots.R: Generate figures of main statistical results
- genPlotsTransmission.R: Generate figures of transmission test results
- getDiagnosisExplained.R: Obtain the ASD diagnosis explained by rare expansions
- largeRepeatTransmissionTest.R: Test of the transmission of large tandem repeats
- mainStats.R: Test for global, gene set and functional element burden of the rare expansions.
- outlierRemovalByEHdnCall.R: Visualize and remove outliers based on the number of EHdn calls per sample
- prepareMapEHdnTRF.R: Map tandem repeats identified by EHdn to known STR from TRF
- testEffectSplicingTSS.R: Test the effect of expansions impacting TSS or splice sites by comparing the loss-of-function intolerance of the genes impacted to other genes.
- TransmittedLociAnalysis.R: Visualize the size of expansions inherited in proband and also transmitted in unaffected siblings.
- XLinkedAnalysis.R: Burden test of x-linked rare expansions.


The following is a list of library files, which contain helper functions used by other scripts in this repository:

- BTlib.py (Python functions)
- BTlib/R/BTlib.R (R functions)
- requiredFunctions.R (R functions)
- gs.asd.2019.RData (R object contains gene sets used in the study)
