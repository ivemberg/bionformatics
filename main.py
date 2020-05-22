import Bio
from Bio import AlignIO
import csv


italy = list(AlignIO.read(open("Results/Kalign/aln-italy.fasta"), "fasta"))

for i in range(1, len(italy)):
    for j in range(len(italy[i].seq)):
       if italy[i].seq[j] != italy[0].seq[j]:
           print("num sequence : " + str(i) + "  position : " + str(j) +
                    " reference base: " + italy[0].seq[j] + " sequence base: " + italy[i].seq[j]) 
            
            

"""           
with open('results-italy.csv', 'w', newline='') as csvfile:
    fields = ['seq name', 'pos', 'ref base', 'seq ref']
    writer = csv.DictWriter(csvfile, fieldnames = fields, delimiter=',')
    writer.writeheader()
    #row = [italy[i].name , j , italy[0].seq[j] , italy[i].seq[j]]
    writer.writerow({'seq name':italy[i].name, 'pos': j, 'ref base': italy[0].seq[j] , 'seq ref': italy[i].seq[j] })
"""
