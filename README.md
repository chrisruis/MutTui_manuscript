# MutTui_manuscript
Data and scripts used in [Calculating and applying pathogen mutational spectra using MutTui](https://www.biorxiv.org/content/10.1101/2023.06.15.545111v1)

The MutTui pipeline can be found [here](https://github.com/chrisruis/MutTui)

### simulations directory
Contains the datasets and scripts used to test whether MutTui reconstructs reliable spectra by simulating 100 datasets from the *Mycobacterium tuberculosis* lineage 4 spectrum and calculating the spectrum. The simulations were run on the lineage 4 tree (MTB_lineage4_rooted.nwk) and [previously calculated](https://www.biorxiv.org/content/10.1101/2022.07.13.499881v1) SBS spectrum (unscaled by genomic composition so the number of mutations is the counts of mutations).

Spectra were simulated through sample_mutations_spectrum.py with "-m 10000" to sample roughly 10,000 mutations. In directory simulated_and_calculated_spectra, the simulated spectrum for each run has suffix "sampled_spectrum.csv" and the spectrum calculated by MutTui has suffix "calculated_spectrum.csv"

### number_of_mutations directory
Contains spectra, converted spectra and downsampled spectra used to assess the number of mutations required to estimated robust spectra for *Mycobacterium tuberculosis* lineage 4, *Streptococcus agalactiae* CC19, SARS-CoV-2 Delta and influenza A H3N2. Subsampling was carried out using subsample_spectra_DNA_RNA.py. Each dataset contains files with suffixes: rescaled.csv (the original SBS spectrum), converted.csv (the SBS spectrum converted to contain 10,000 mutations in the same proportions as the original spectrum in rescaled.csv), cosine.csv (the number of mutations and cosine similarity between each downsampled spectrum and the original spectrum)

### tree_topology directory
Contains datasets and MutTui outputs used to assess impact of tree topology on spectrum calculation. In each case, MutTui was run on 10 tree topologies sampled from a posterior distribution. The S_agalactiae_CC1 and MERS directories each contain the alignment, reference and position conversion file (same for each run) and 10 run directories that contain the respective tree and output from MutTui

### tree_rooting directory
Contains datasets and MutTui output used to assess impact of rooting strategy on the *Streptococcus agalactiae* CC1 spectrum. The alignment, reference and position conversion file were the same for each run. The outgroup_rooting, midpoint_rooting and date_rooting directories contain the phylogenetic tree and MutTui output for each rooting strategy


### nth_gene_signatures directory
Contains the files and spectra used to calculate the *nth* gene signature in *Mycobacterium leprae* and *Mycobacterium abscessus*.

**Mycobacterium leprae**
The *M. leprae nth* gene signature was calculated from a phylogenetic tree that includes several hypermutator lineages. The mycobacterium_leprae directory includes the alignment (M_leprae.fasta), unlabelled tree (M_leprae.nwk), position conversion file (conversion.txt) and reference (reference.fasta). The tree was labelled to identify the hypermutator lineages using MutTui. Initially, node labels were added to the unlabelled tree to identify the branches on which the label changes using:
```
MutTui label-nodes -t M_leprae.nwk -o M_leprae_nodes_labelled.nwk
```

By viewing M_leprae_nodes_labelled.nwk with node labels (for example using FigTree), there are 4 tip branches with evidence of hypermutation and one internal branch leading to two sequences. The node labels show the internal branch is Node123). Therefore the tree can be labelled to separate the hypermutators (label nth) from the background branches (label B) using:
```
MutTui label-tree -t M_leprae.nwk -r B -s SRR6241727____nth SRR6241736_1____nth SRR6241800____nth SRR6241758____nth Node123____nth -o M_leprae_nodes_labelled.nwk
```

Then MutTui can be run to calculate the spectrum of the hypermutator lineages (and background branches separately) using:
```
MutTui run -a M_leprae.fasta -t M_leprae.nwk -lt M_leprae_labelled.nwk -c conversion.txt -r reference.fasta -o muttui_out --include_all_branches
```

The *nth* gene signature is in mutational_spectrum_label_nth_rescaled.csv

**Mycobacterium abscessus**
The *M. abscessus nth* gene signature was calculated using SNPs from deep sequencing data of a patient infection where *nth* was knocked out. The mutations are in all_mutations.csv and the spectrum of these mutations can be calculated using:
```
MutTui korimuto -v all_mutations.csv -r reference.fa --variant -o nth
```

### environment_lung_comparison directory
Contains data used to compare environmental and lung spectra in *Mycobacteria* and *Burkholderia*. The alignments (named with the clade name followed by .fasta), phylogenetic trees (named with the clade name followed by .nwk), position conversion files (named conversion.txt), references (named reference.fasta) and SBS spectra (named with the clade name followed by \_SBS.csv) for each dataset are in the respective directory. *M. canettii* and *M. kansasii* also contain labelled trees (named with the clade name followed by \_labelled.nwk) as their SBS spectra are calculated on a subset of branches.

Directory spectrum_comparisons contains spectrum clustering and mutation type comparisons. SBS_combined.csv contains all included SBS spectra. Coordinates of PCA based on overall spectrum compositions are in SBS_PCA_coordinates.csv. The proportions of C>A mutations and T>C mutations in each spectrum are in CA_vs_TC.csv
