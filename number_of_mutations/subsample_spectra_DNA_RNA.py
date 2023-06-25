#Subsampling analyses with cosine similarity comparisons for 2 DNA and 2 RNA datasets

from scipy import spatial
import random

#Sums mutation types from a spectrum
def sumMT(spectrum, rna=False):
    mt = list()

    mt.append(sum(spectrum[:16]))
    mt.append(sum(spectrum[16:32]))
    mt.append(sum(spectrum[32:48]))
    mt.append(sum(spectrum[48:64]))
    mt.append(sum(spectrum[64:80]))

    if rna:
        mt.append(sum(spectrum[80:96]))
        mt.append(sum(spectrum[96:112]))
        mt.append(sum(spectrum[112:128]))
        mt.append(sum(spectrum[128:144]))
        mt.append(sum(spectrum[144:160]))
        mt.append(sum(spectrum[160:176]))
        mt.append(sum(spectrum[176:]))
    else:
        mt.append(sum(spectrum[80:]))

    return(mt)

#Count mutations in input spectrum
def mutationCounts(spectrum, rna=False):
    #Set of mutations
    if rna:
        mutations = ["A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T","A[G>T]A","A[G>T]C","A[G>T]G","A[G>T]T","C[G>T]A","C[G>T]C","C[G>T]G","C[G>T]T","G[G>T]A","G[G>T]C","G[G>T]G","G[G>T]T","T[G>T]A","T[G>T]C","T[G>T]G","T[G>T]T","A[G>C]A","A[G>C]C","A[G>C]G","A[G>C]T","C[G>C]A","C[G>C]C","C[G>C]G","C[G>C]T","G[G>C]A","G[G>C]C","G[G>C]G","G[G>C]T","T[G>C]A","T[G>C]C","T[G>C]G","T[G>C]T","A[G>A]A","A[G>A]C","A[G>A]G","A[G>A]T","C[G>A]A","C[G>A]C","C[G>A]G","C[G>A]T","G[G>A]A","G[G>A]C","G[G>A]G","G[G>A]T","T[G>A]A","T[G>A]C","T[G>A]G","T[G>A]T","A[A>T]A","A[A>T]C","A[A>T]G","A[A>T]T","C[A>T]A","C[A>T]C","C[A>T]G","C[A>T]T","G[A>T]A","G[A>T]C","G[A>T]G","G[A>T]T","T[A>T]A","T[A>T]C","T[A>T]G","T[A>T]T","A[A>G]A","A[A>G]C","A[A>G]G","A[A>G]T","C[A>G]A","C[A>G]C","C[A>G]G","C[A>G]T","G[A>G]A","G[A>G]C","G[A>G]G","G[A>G]T","T[A>G]A","T[A>G]C","T[A>G]G","T[A>G]T","A[A>C]A","A[A>C]C","A[A>C]G","A[A>C]T","C[A>C]A","C[A>C]C","C[A>C]G","C[A>C]T","G[A>C]A","G[A>C]C","G[A>C]G","G[A>C]T","T[A>C]A","T[A>C]C","T[A>C]G","T[A>C]T"]
    else:
        mutations = ["A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"]

    s = list()
    with open(spectrum) as f:
        next(f)
        for l in f:
            for i in [l.strip().split(",")[0]] * int(l.strip().split(",")[1]):
                s.append(i)
    #Mutation counts in original spectrum
    sC = list()
    for m in mutations:
        sC.append(s.count(m))
    
    #Extract mutation types
    mt = sumMT(sC, rna)
    
    return(s, sC, mt)

#Randomly subsample the mutation set at different levels, calculate the spectrum and compare with the original
def subsampleSpectrum(s, c, m, o, rna=False):
    #Set of mutations
    if rna:
        mutations = ["A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T","A[G>T]A","A[G>T]C","A[G>T]G","A[G>T]T","C[G>T]A","C[G>T]C","C[G>T]G","C[G>T]T","G[G>T]A","G[G>T]C","G[G>T]G","G[G>T]T","T[G>T]A","T[G>T]C","T[G>T]G","T[G>T]T","A[G>C]A","A[G>C]C","A[G>C]G","A[G>C]T","C[G>C]A","C[G>C]C","C[G>C]G","C[G>C]T","G[G>C]A","G[G>C]C","G[G>C]G","G[G>C]T","T[G>C]A","T[G>C]C","T[G>C]G","T[G>C]T","A[G>A]A","A[G>A]C","A[G>A]G","A[G>A]T","C[G>A]A","C[G>A]C","C[G>A]G","C[G>A]T","G[G>A]A","G[G>A]C","G[G>A]G","G[G>A]T","T[G>A]A","T[G>A]C","T[G>A]G","T[G>A]T","A[A>T]A","A[A>T]C","A[A>T]G","A[A>T]T","C[A>T]A","C[A>T]C","C[A>T]G","C[A>T]T","G[A>T]A","G[A>T]C","G[A>T]G","G[A>T]T","T[A>T]A","T[A>T]C","T[A>T]G","T[A>T]T","A[A>G]A","A[A>G]C","A[A>G]G","A[A>G]T","C[A>G]A","C[A>G]C","C[A>G]G","C[A>G]T","G[A>G]A","G[A>G]C","G[A>G]G","G[A>G]T","T[A>G]A","T[A>G]C","T[A>G]G","T[A>G]T","A[A>C]A","A[A>C]C","A[A>C]G","A[A>C]T","C[A>C]A","C[A>C]C","C[A>C]G","C[A>C]T","G[A>C]A","G[A>C]C","G[A>C]G","G[A>C]T","T[A>C]A","T[A>C]C","T[A>C]G","T[A>C]T"]
    else:
        mutations = ["A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"]

    #Levels of subsampling
    levels = [10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 50]

    #Lists of cosine similarities at each level
    sim = dict()
    mtSim = dict()
    
    for l in levels:
        print("Level:", l)
        sim[l] = list()
        mtSim[l] = list()

        for i in range(1000):
            #Subsample mutations
            ss = random.choices(s, k = l)

            sSpectrum = list()
            #Count each mutation in subsample
            for eM in mutations:
                sSpectrum.append(ss.count(eM))
            #Sum mutation types
            sMT = sumMT(sSpectrum, rna)
        
            #Calculate distance to original spectrum
            sim[l].append(float(1) - spatial.distance.cosine(c, sSpectrum))
            mtSim[l].append(float(1) - spatial.distance.cosine(m, sMT))
    
    #Write similarities
    outFile = open(o + ".csv", "w")
    outFile.write("Level,Cosine_similarity\n")
    outMT = open(o + "_mt.csv", "w")
    outMT.write("Level,Cosine_similarity\n")
    for l in levels:
        for eS in sim[l]:
            outFile.write(str(l) + "," + str(eS) + "\n")
        for eS in mtSim[l]:
            outMT.write(str(l) + "," + str(eS) + "\n")
    outFile.close()

#Import spectra and calculate mutation counts
mtb, mtbC, mtbMT = mutationCounts("MTB_L4_converted.csv")
cc19, cc19C, cc19MT = mutationCounts("CC19_converted.csv")
delta, deltaC, deltaMT = mutationCounts("delta_converted.csv", rna=True)
h3n2, h3n2C, h3n2MT = mutationCounts("H3N2_converted.csv", rna=True)

#Iterate through the spectra, subsample and compare with original
subsampleSpectrum(mtb, mtbC, mtbMT, "MTB_L4_cosine")
subsampleSpectrum(cc19, cc19C, cc19MT, "CC19_cosine")
subsampleSpectrum(delta, deltaC, deltaMT, "delta_cosine", rna=True)
subsampleSpectrum(h3n2, h3n2C, h3n2MT, "H3N2_cosine", rna=True)