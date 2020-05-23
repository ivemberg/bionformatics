import Bio
from Bio import AlignIO
import csv
import glob

# Full version

for f in glob.glob('Results/Kalign/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        print("sequence : " + aln[i].id, file=open("Seq_differences/full/kalign-output.txt", "a"))
        for j in range(len(aln[i].seq)):
            if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                print("  position : " + str(j) + " reference base: " + aln[0].seq[j] + " sequence base: " + aln[i].seq[j], file=open("Seq_differences/full/kalign-output.txt", "a")) 

for f in glob.glob('Results/Clustal Omega/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        print("sequence : " + aln[i].id, file=open("Seq_differences/full/clustal-output.txt", "a"))
        for j in range(len(aln[i].seq)):
            if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                print("  position : " + str(j) + " reference base: " + aln[0].seq[j] + " sequence base: " + aln[i].seq[j], file=open("Seq_differences/full/clustal-output.txt", "a")) 

for f in glob.glob('Results/MAFFT/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        print("sequence : " + aln[i].id, file=open("Seq_differences/full/MAFFT-output.txt", "a"))
        for j in range(len(aln[i].seq)):
            if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                print("  position : " + str(j) + " reference base: " + aln[0].seq[j] + " sequence base: " + aln[i].seq[j], file=open("Seq_differences/full/MAFFT-output.txt", "a")) 

haln = list(AlignIO.read("Results/Horizontal/aln-clustal-horizontal.fasta", "fasta"))
for i in range(1, len(haln)):
    print("sequence : " + haln[i].id, file=open("Seq_differences/full/clustal-horizontal-output.txt", "a"))
    for j in range(len(haln[i].seq)):
        if haln[i].seq[j].upper() != haln[0].seq[j].upper():
            print("  position : " + str(j) + " reference base: " + haln[0].seq[j] + " sequence base: " + haln[i].seq[j], file=open("Seq_differences/full/clustal-horizontal-output.txt", "a")) 

haln = list(AlignIO.read("Results/Horizontal/aln-kalign-horizontal.fasta", "fasta"))
for i in range(1, len(haln)):
    print("sequence : " + haln[i].id, file=open("Seq_differences/full/kalign-horizontal-output.txt", "a"))
    for j in range(len(haln[i].seq)):
        if haln[i].seq[j].upper() != haln[0].seq[j].upper():
            print("  position : " + str(j) + " reference base: " + haln[0].seq[j] + " sequence base: " + haln[i].seq[j], file=open("Seq_differences/full/kalign-horizontal-output.txt", "a")) 

haln = list(AlignIO.read("Results/Horizontal/aln-MAFFT-horizontal.fasta", "fasta"))
for i in range(1, len(haln)):
    print("sequence : " + haln[i].id, file=open("Seq_differences/full/MAFFT-horizontal-output.txt", "a"))
    for j in range(len(haln[i].seq)):
        if haln[i].seq[j].upper() != haln[0].seq[j].upper():
            print("  position : " + str(j) + " reference base: " + haln[0].seq[j] + " sequence base: " + haln[i].seq[j], file=open("Seq_differences/full/MAFFT-horizontal-output.txt", "a")) 