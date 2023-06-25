# MutTui_manuscript
Data and scripts used in [Calculating and applying pathogen mutational spectra using MutTui](https://www.biorxiv.org/content/10.1101/2023.06.15.545111v1)

### simulations directory
Contains the datasets and scripts used to test whether MutTui reconstructs reliable spectra by simulating 100 datasets from the Mycobacterium tuberculosis lineage 4 spectrum and calculating the spectrum. The simulations were run on the lineage 4 tree (MTB_lineage4_rooted.nwk) and [previously calculated](https://www.biorxiv.org/content/10.1101/2022.07.13.499881v1) SBS spectrum (unscaled by genomic composition so the number of mutations is the counts of mutations).

Spectra were simulated through sample_mutations_spectrum.py with "-m 10000" to sample roughly 10,000 mutations. In directory simulated_and_calculated_spectra, the simulated spectrum for each run has suffix "sampled_spectrum.csv" and the spectrum calculated by MutTui has suffix "calculated_spectrum.csv"

### number_of_mutations directory
Contains spectra, converted spectra and downsampled spectra used to assess the number of mutations required to estimated robust spectra for Mycobacterium tuberculosis lineage 4, Streptococcus agalactiae CC19, SARS-CoV-2 Delta and influenza A H3N2. Subsampling was carried out using subsample_spectra_DNA_RNA.py. Each dataset contains files with suffixes: rescaled.csv (the original SBS spectrum), converted.csv (the SBS spectrum converted to contain 10,000 mutations in the same proportions as the original spectrum in rescaled.csv), cosine.csv (the number of mutations and cosine similarity between each downsampled spectrum and the original spectrum)

### tree_topology directory
Contains datasets and MutTui outputs used to assess impact of tree topology on spectrum calculation. In each case, MutTui was run on 10 tree topologies sampled from a posterior distribution. The S_agalactiae_CC1 and MERS directories each contain the alignment, reference and position conversion file (same for each run) and 10 run directories that contain the respective tree and output from MutTui

### tree_rooting directory
Contains datasets and MutTui output used to assess impact of rooting strategy on the Streptococcus agalactiae CC1 spectrum. The alignment, reference and position conversion file were the same for each run. The outgroup_rooting, midpoint_rooting and date_rooting directories contain the phylogenetic tree and MutTui output for each rooting strategy
