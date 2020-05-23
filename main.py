import Bio
from Bio import AlignIO
import csv
import glob

fasta_list = glob()
italy = list(AlignIO.read(open("Results/Clustal Omega/italy.fa"), "fasta"))

for i in range(1, len(italy)):
    for j in range(len(italy[i].seq)):
       if italy[i].seq[j] != italy[0].seq[j]:
                print("ITALY \n" + "sequence : " + italy[i].id.partition("|")[2] + "  position : " + str(j) +
                    " reference base: " + italy[0].seq[j] + " sequence base: " + italy[i].seq[j], file=open("output.txt", "a")) 

"""               
ny = list(AlignIO.read(open("Results/Kalign/aln-newyork.fasta"), "fasta"))

for i in range(1, len(ny)):
    for j in range(len(ny[i].seq)):
       if ny[i].seq[j] != ny[0].seq[j]:
                print("NY \n" + "sequence : " + ny[i].id.partition("|")[2] + "  position : " + str(j) +
                    " reference base: " + ny[0].seq[j] + " sequence base: " + ny[i].seq[j], file=open("output.txt", "a")) 

russia = list(AlignIO.read(open("Results/Kalign/aln-russia.fasta"), "fasta"))

for i in range(1, len(russia)):
    for j in range(len(russia[i].seq)):
       if russia[i].seq[j] != russia[0].seq[j]:
                print("RUSSIA \n" + "sequence : " + russia[i].id.partition("|")[2] + "  position : " + str(j) +
                    " reference base: " + russia[0].seq[j] + " sequence base: " + russia[i].seq[j], file=open("output.txt", "a")) 

spain = list(AlignIO.read(open("Results/Kalign/aln-spain.fasta"), "fasta"))

for i in range(1, len(spain)):
    for j in range(len(spain[i].seq)):
       if spain[i].seq[j] != spain[0].seq[j]:
                print("SPAIN \n" + "sequence : " + spain[i].id.partition("|")[2] + "  position : " + str(j) +
                    " reference base: " + spain[0].seq[j] + " sequence base: " + spain[i].seq[j], file=open("output.txt", "a")) 

"""

"""           
with open('results-italy.csv', 'w+', newline='') as csvfile:
    fields = ['seq name', 'pos', 'ref base', 'seq ref']
    writer = csv.DictWriter(csvfile, fieldnames = fields, delimiter=',')
    writer.writeheader()
    #row = [italy[i].name , j , italy[0].seq[j] , italy[i].seq[j]]
    writer.writerow({'seq name':italy[i].name, 'pos': j, 'ref base': italy[0].seq[j] , 'seq ref': italy[i].seq[j] })
"""
