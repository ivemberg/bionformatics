import Bio
from Bio import GenBank, SeqIO, AlignIO
import glob

"""
t1 = find_genes(b-1, ref, seq)
                        t2 = find_genes(b-2, ref, seq)
                        if t1:
                            position.append(t1[1])
                        if t2:
                            position.append(t2[1])
"""

def find_genes(b,ref,seq):
    result = []
    for gb_record in SeqIO.parse(open("sequence.gb","r"), "genbank") :
        for features in gb_record.features:
            if features.type == "CDS":
                for i in range(min(features.location), max(features.location) + 1, 3): #265 268 271 (inizio dei codoni)
                    cod = ref[i]+ref[i+1]+ref[i+2]
                    m_cod = seq[i]+seq[i+1]+seq[i+2]
                    if b in range(i, i+3): #se b è 79 80 81
                        element = [b, cod, m_cod.upper()]
                        result.append(element)
                 
    if result:
        for i in range(len(result)):
            for j in range(i+1, len(result)):
                if (result[i][0] == (result[j][0]- 1)) and (result[i][1] == result[j][1]):
                    result.remove(result[i+1])
        print(result)


                    


                        
    """
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

    """


def main():
    
    """"
    for f in glob.glob('Results/Kalign/aln-russia.fasta'):
        aln = list(AlignIO.read(open(f), "fasta"))
        for i in range(1, len(aln)):
            for j in range(len(aln[i].seq)):
                if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                    find_genes(j, aln[0], aln[i])
    
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
                        s = ""
                        r = ""
                        l = 0
                else:
                    if s != "" and r != "":
                        element =  [aln[i].id.split('/')[0], j - l, l, r, s]
                        diffs = []
                        k.append(element)
                        for position in range(j - l, j+1) :
                            if aln[i].seq[position].upper() != aln[0].seq[position].upper():
                                diff = find_genes(position, aln[0], aln[i])
                                if diff:
                                    diffs.append(find_genes(position, aln[0], aln[i]))
                        #print(element)
                        if diffs:
                            print(diffs)
                        print("")
                        s = ""
                        r = ""
                        l = 0
"""
    s = ""
    r = ""
    k = []
    l=0
    for f in glob.glob('Results/Kalign/aln-russia.fasta'):
        aln = list(AlignIO.read(open(f), "fasta"))
        for i in range(1, len(aln)):
            for j in range(len(aln[i].seq)):
                if aln[i].seq[j].upper() != aln[0].seq[j].upper():
                    find_genes(j,aln[0],aln[i])
                    s += aln[i].seq[j]
                    r += aln[0].seq[j]
                    l += 1
                    if (j == len(aln[i].seq)-1):
                        element =  [aln[i].id.split('/')[0], j - l, l, r, s]
                        k.append(element)
                        s = ""
                        r = ""
                        l = 0
                else:
                    if s != "" and r != "":
                        element =  [aln[i].id.split('/')[0], j - l, l, r, s]
                        k.append(element)
                        s = ""
                        r = ""
                        l = 0


if __name__ == "__main__":
	main()



