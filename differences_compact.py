import Bio
from Bio import AlignIO
import csv
import glob

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
            else:
                if s != "" and r != "":
                    element =  [aln[i].id, j - l, r, s]
                    k.append(element)
                    s = ""
                    r = ""
                    l = 0

with open('Seq_differences/full/kalign-compact-output.txt', 'w') as f:
    for item in k:
        f.write("%s\n" % item)        

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
            else:
                if s != "" and r != "":
                    element =  [aln[i].id, j - l, r, s]
                    c.append(element)
                    s = ""
                    r = ""
                    l = 0

with open('Seq_differences/full/clustal-omega-compact-output.txt', 'w') as f:
    for item in c:
        f.write("%s\n" % item)       

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
            else:
                if s != "" and r != "":
                    element =  [aln[i].id, j - l, r, s]
                    m.append(element)
                    s = ""
                    r = ""
                    l = 0

with open('Seq_differences/full/MAFFT-compact-output.txt', 'w') as f:
    for item in m:
        f.write("%s\n" % item)    

s = ""
r = ""
l = 0
hk = []
haln = list(AlignIO.read("Results/Horizontal/aln-kalign-horizontal.fasta", "fasta"))
for i in range(1, len(haln)):
        for j in range(len(haln[i].seq)):
            if haln[i].seq[j].upper() != haln[0].seq[j].upper():
                s += haln[i].seq[j]
                r += haln[0].seq[j]
                l += 1
            else:
                if s != "" and r != "":
                    element =  [haln[i].id, j - l, r, s]
                    hk.append(element)
                    s = ""
                    r = ""
                    l = 0

with open('Seq_differences/full/kalign-compact-horizontal-output.txt', 'w') as f:
    for item in hk:
        f.write("%s\n" % item)   

s = ""
r = ""
l = 0
hc = []
haln = list(AlignIO.read("Results/Horizontal/aln-clustal-horizontal.fasta", "fasta"))
for i in range(1, len(haln)):
        for j in range(len(haln[i].seq)):
            if haln[i].seq[j].upper() != haln[0].seq[j].upper():
                s += haln[i].seq[j]
                r += haln[0].seq[j]
                l += 1
            else:
                if s != "" and r != "":
                    element =  [haln[i].id, j - l, r, s]
                    hc.append(element)
                    s = ""
                    r = ""
                    l = 0

with open('Seq_differences/full/clustal-compact-horizontal-output.txt', 'w') as f:
    for item in hc:
        f.write("%s\n" % item)  

s = ""
r = ""
l = 0
hm = []
haln = list(AlignIO.read("Results/Horizontal/aln-MAFFT-horizontal.fasta", "fasta"))
for i in range(1, len(haln)):
        for j in range(len(haln[i].seq)):
            if haln[i].seq[j].upper() != haln[0].seq[j].upper():
                s += haln[i].seq[j]
                r += haln[0].seq[j]
                l += 1
            else:
                if s != "" and r != "":
                    element =  [haln[i].id, j - l, r, s]
                    hm.append(element)
                    s = ""
                    r = ""
                    l = 0

with open('Seq_differences/full/MAFFT-compact-horizontal-output.txt', 'w') as f:
    for item in hm:
        f.write("%s\n" % item)  