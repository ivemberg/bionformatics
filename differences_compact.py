import Bio
from Bio import AlignIO
import csv
import glob

s = ""
r = ""
c = []
l = 0

for f in glob.glob('Results/Kalign/*.fasta'):
    aln = list(AlignIO.read(open(f), "fasta"))
    for i in range(1, len(aln)):
        # print("sequence : " + aln[i].id, file=open("Seq_differences/full/kalign-output.txt", "a"))
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

with open('Seq_differences/full/kalign-output.txt', 'w') as f:
    for item in c:
        f.write("%s\n" % item)        


