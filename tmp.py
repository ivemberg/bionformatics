import Bio
from Bio import AlignIO
import csv
import glob

s = ""
r = ""
l = 0
k = []

aln = list(AlignIO.read("Results/Kalign/aln-italy.fasta", "fasta"))
for i in range(5,6): #considero allineamento 5
    for j in range(len(aln[i].seq)): 
        if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                s += aln[i].seq[j]
                r += aln[0].seq[j]
                l += 1
        else:
            if s != "" and r != "":
                element =  [aln[i].id, j - l, r, s]
                k.append(element)
                s = ""
                r = ""
                l = 0
    

with open('tmp-kalign-c.txt', 'w') as f:
    for item in k:
        f.write("%s\n" % item)      
