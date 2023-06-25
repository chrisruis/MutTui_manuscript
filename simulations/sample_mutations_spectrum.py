#Simulates mutations on a phylogenetic tree based on a given spectrum
#Samples the number of mutations on each branch from a Poisson distribution with lambda equal to branch length
#Randomly samples that number of mutations from the spectrum

from Bio import Phylo
from numpy.random import poisson
from random import sample
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", help = "Newick tree")
    parser.add_argument("-s", help = "Mutational spectrum")
    parser.add_argument("-m", help = "Target number of mutations, branch lengths will be scaled to result in this number of expected mutations, default = 10000", default = "10000")
    parser.add_argument("-o", help = "Prefix of output files. 3 files are written: an alignment of all sequences named with this prefix followed by .fasta, a reference genome for use with MutTui named with this prefix followed by _reference.fasta, the sampled spectrum of mutations named with this prefix followed by _sampled_spectrum.csv")
    args = parser.parse_args()

    #96 mutations
    sMutations = ["A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"]

    #Import tree
    tree = Phylo.read(args.t, "newick")

    #Import spectrum and create empty spectrum
    s = list()
    spectrum = dict()
    with open(args.s) as f:
        next(f)
        for l in f:
            s += [l.strip().split(",")[0]] * int(l.strip().split(",")[1])
            spectrum[l.strip().split(",")[0]] = 0

    #Calculate total branch length in tree
    tbl = float(0)
    for c in tree.find_clades():
        #Don't analyse the root
        if c.branch_length:
            tbl += c.branch_length
    #Calculate scale factor for branches
    sf = float(args.m)/tbl

    #Mutations with positions and sequences
    mutations = dict()

    #Iterator for mutations
    i = 1

    allM = 0

    #Iterate through the branches and sample mutations
    for c in tree.find_clades():
        #Don't analyse the root or the branches immediately downstream as these are excluded by MutTui
        if len(tree.get_path(c)) >= 2:
            #Sample mutations on branch
            nm = poisson(c.branch_length * sf)
            allM += nm
            for eM in range(nm):
                sM = sample(s, 1)[0]
                spectrum[sM] += 1
                mutations[i] = list()
                mutations[i].append(sM)
                mutations[i].append(set())
                for tip in c.get_terminals():
                    mutations[i][1].add(tip.name)
                i += 1
    
    #Sampled spectrum
    ss = dict()
    for eM in sMutations:
        ss[eM] = 0
    
    #Combine sequences
    sequences = dict()
    for tip in tree.get_terminals():
        sequences[tip.name] = ""
    for m in mutations:
        mut = mutations[m][0]
        ss[mut] += 1
        for tip in tree.get_terminals():
            if tip.name in mutations[m][1]:
                sequences[tip.name] += mut[0] + mut[4] + mut[6]
            else:
                sequences[tip.name] += mut[0] + mut[2] + mut[6]

    #Used to write reference
    firstSeq = True

    #Write sequences
    out = open(args.o + ".fasta", "w")
    outRef = open(args.o + "_reference.fasta", "w")
    for eS in sequences:
        if firstSeq:
            outRef.write(">" + eS + "\n" + sequences[eS] + "\n")
            firstSeq = False
        out.write(">" + eS + "\n" + sequences[eS] + "\n")
    out.close()
    outRef.close()

    #Write spectrum
    outS = open(args.o + "_sampled_spectrum.csv", "w")
    outS.write("Substitution,Number_of_mutations\n")
    for eMut in ss:
        outS.write(eMut + "," + str(ss[eMut]) + "\n")
    outS.close()