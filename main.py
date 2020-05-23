import Bio
from Bio import AlignIO
import csv
import glob


for f in glob.glob('Results/Kalign/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        print("sequence : " + aln[i].id, file=open("kalign-output.txt", "a"))
        for j in range(len(aln[i].seq)):
            if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                print("  position : " + str(j) + " reference base: " + aln[0].seq[j] + " sequence base: " + aln[i].seq[j], file=open("kalign-output.txt", "a")) 

"""
for f in glob.glob('Results/Clustal Omega/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        print("sequence : " + aln[i].id, file=open("clustalomega-output.txt", "a"))
        for j in range(len(aln[i].seq)):
            if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                print("  position : " + str(j) + " reference base: " + aln[0].seq[j] + " sequence base: " + aln[i].seq[j], file=open("clustalomega-output.txt", "a")) 
"""

"""
for f in glob.glob('Results/MAFFT/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        print("sequence : " + aln[i].id, file=open("mafft-output.txt", "a"))
        for j in range(len(aln[i].seq)):
            if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                print("  position : " + str(j) + " reference base: " + aln[0].seq[j] + " sequence base: " + aln[i].seq[j], file=open("mafft-output.txt", "a")) 
"""

