import Bio
from Bio import GenBank, SeqIO, AlignIO
import glob

def find_genes(b,ref,seq):
    result = []
    gb_file = "sequence.gb"
    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        for features in gb_record.features :
            if (features.type == "CDS"): 
                if b in features.location: # intervallo CDS
                    p = b - min(features.location) + 1 
                    if p%3 == 0:
                        cod = ref[b-2]+ref[b-1]+ref[b] #codone originale
                        m_cod = seq[b-2]+seq[b-1]+seq[b]
                    elif p%3 == 1 and (b+2) < max(features.location):
                        cod = ref[b]+ref[b+1]+ref[b+2]
                        m_cod = seq[b]+seq[b+1]+seq[b+2]
                    elif p%3 == 2 and (b+1) < max(features.location):
                        cod = ref[b-1]+ref[b]+ref[b+1]
                        m_cod = seq[b-1]+seq[b]+seq[b+1]
                    e = seq.id, b, features.qualifiers['db_xref'][0], cod, m_cod, min(features.location), max(features.location)
                    result.append(e)
    if result:
        print(result)


def main():
    
    for f in glob.glob('Results/Kalign/aln-russia.fasta'):
        aln = list(AlignIO.read(open(f), "fasta"))
        for i in range(1, len(aln)):
            for j in range(len(aln[i].seq)):
                if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                    find_genes(j, aln[0], aln[i])
    """
    s = ""
    r = ""
    l = 0
    k = []
    for f in glob.glob('Results/Kalign/aln-russia.fasta'):
        aln = list(AlignIO.read(open(f), "fasta"))
        for i in range(1, len(aln)):
            for j in range(len(aln[i].seq)):
                if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                    s += aln[i].seq[j]
                    r += aln[0].seq[j]
                    l += 1
                    if (j == len(aln[i].seq)-1):
                        element =  [aln[i].id.split('/')[0], j - l, l, r, s]
                        k.append(element)
                        find_genes(j-l, aln[0], aln[i])
                        s = ""
                        r = ""
                        l = 0
                else:
                    if s != "" and r != "":
                        element =  [aln[i].id.split('/')[0], j - l, l, r, s]
                        k.append(element)
                        find_genes(j-l, aln[0], aln[i])
                        s = ""
                        r = ""
                        l = 0
"""

if __name__ == "__main__":
	main()



