# MutTui_manuscript
Data and scripts used in [Calculating and applying pathogen mutational spectra using MutTui](https://www.biorxiv.org/content/10.1101/2023.06.15.545111v1)

### simulations directory
Contains the datasets and scripts used to test whether MutTui reconstructs reliable spectra by simulating 100 datasets from the Mycobacterium tuberculosis lineage 4 spectrum and calculating the spectrum. The simulations were run on the lineage 4 tree (MTB_lineage4_rooted.nwk) and [previously calculated](https://www.biorxiv.org/content/10.1101/2022.07.13.499881v1) SBS spectrum (unscaled by genomic composition so the number of mutations is the counts of mutations).

Spectra were simulated through sample_mutations_spectrum.py with "-m 10000" to sample roughly 10,000 mutations. In directory simulated_and_calculated_spectra, the simulated spectrum for each run has suffix "sampled_spectrum.csv" and the spectrum calculated by MutTui has suffix "calculated_spectrum.csv"
