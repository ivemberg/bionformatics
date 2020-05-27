import Bio
from Bio import AlignIO
import csv
import glob
import collections

s = ""
r = ""
l = 0

k = []
for f in glob.glob('Results/Kalign/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        for j in range(len(aln[i].seq)):
            if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                s += aln[i].seq[j]
                r += aln[0].seq[j]
                l += 1
                if (j == len(aln[i].seq)-1):
                    element =  aln[i].id.split('/')[0], j - l, l, r, s
                    k.append(element)
                    s = ""
                    r = ""
                    l = 0
            else:
                if s != "" and r != "":
                    element =  aln[i].id.split('/')[0], j - l, l, r, s
                    k.append(element)
                    s = ""
                    r = ""
                    l = 0

s = ""
r = ""
l = 0
c = []
for f in glob.glob('Results/Clustal Omega/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        for j in range(len(aln[i].seq)):
            if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                s += aln[i].seq[j]
                r += aln[0].seq[j]
                l += 1
                if (j == len(aln[i].seq)-1):
                    element =  aln[i].id.split('/')[0], j - l, l, r, s
                    c.append(element)
                    s = ""
                    r = ""
                    l = 0
            else:
                if s != "" and r != "":
                    element =  aln[i].id.split('/')[0], j - l, l, r, s
                    c.append(element)
                    s = ""
                    r = ""
                    l = 0

s = ""
r = ""
l = 0
m = []
for f in glob.glob('Results/MAFFT/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        for j in range(len(aln[i].seq)):
            if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                s += aln[i].seq[j]
                r += aln[0].seq[j]
                l += 1
                if (j == len(aln[i].seq)-1):
                    element =  aln[i].id.split('/')[0], j - l, l, r.upper(), s.upper()
                    m.append(element)
                    s = ""
                    r = ""
                    l = 0
            else:
                if s != "" and r != "":
                    element =  aln[i].id.split('/')[0], j - l, l, r.upper(), s.upper()
                    m.append(element)
                    s = ""
                    r = ""
                    l = 0

"""
print(len(k), len(c), len(m))
for i in range(len(k)):
    if (k[i] != c[i]):
        print("Kalign: ", k[i], "Clustal Omega: ", c[i], "\n", file = open("tool-diff.txt", "a"))
    if (m[i] != c[i]):
        print("MAFFT: ", k[i], "Clustal Omega: ", c[i], "\n", file = open("tool-diff.txt", "a"))
    if (m[i] != k[i]):
        print("MAFFT: ", k[i], "Kalign: ", c[i], "\n", file = open("tool-diff.txt", "a"))                
"""

# k, c e m sono list
# print(k[0][0]) id sequenza k[0][0][2:] senza 0_
# print(k[0][1]) pos
# print(k[0][2]) length
# print(k[0][3]) ref
# print(k[0][4]) mut

"""
non funziona 

k_clean = []
it = []
for i in range(len(k)):
    if ("Italy" in k[i][0]):
        it = k[0][1], k[0][2], k[0][3], k[0][4]
        k_clean.append(it)
print(k_clean)
print([item for item, count in collections.Counter(k_clean).items() if count > 1], file = open("kalign-duplicates.txt", "a"))



a = [1,2,3,2,1,5,6,5,5,5]

import collections
print([item for item, count in collections.Counter(a).items() if count > 1])

## [1, 2, 5]
"""
     