import Bio
from Bio import AlignIO
import csv
import glob

# check differences in tool output 
# check similarities in mutations (for the same country)
# check similarities in horizontal alignment

# create list instead of []

# vertical 
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


tool_d = []

if len(c) == len(k) == len(m):
    for i in range(len(k)):
        if not(k[i] == c[i] == m[i]): # elemento della lista non uguale
            d = "K: ", k[i], "C: ", c[i], "M: ", m[i]
            tool_d.append(d)

print(tool_d, file = open("tool-diff.txt", "w")) # non sono sicura funzioni 


