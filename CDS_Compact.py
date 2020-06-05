import Bio
from Bio import GenBank, SeqIO, AlignIO
import glob

def translate(codon): 
      
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    return table[codon]

def find_genes(b,ref,seq):
    gb_file = "sequence.gb"
    result = []
    ref = ref.upper()
    seq = seq.upper()

    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        for features in gb_record.features :
            if (features.type == "CDS"): 
                if b in features.location: # intervallo CDS
                    p = b - min(features.location) + 1 
                    positions = []
                    
                    if p%3 == 0:
                        cod = ref[b-2]+ref[b-1]+ref[b] #codone originale
                        m_cod = seq[b-2]+seq[b-1]+seq[b]
                        # Check adiacent positions
                        if ref[b-2] != seq[b-2]:
                            positions.append(b-2)
                        if ref[b-1] != seq[b-1]:
                            positions.append(b-1)
                        positions.append(b)

                    elif p%3 == 1 and (b+2) < max(features.location):
                        cod = ref[b]+ref[b+1]+ref[b+2]
                        m_cod = seq[b]+seq[b+1]+seq[b+2]
                        positions.append(b)
                        if ref[b+1] != seq[b+1]:
                            positions.append(b+1)
                        if ref[b+2] != seq[b+2]:
                            positions.append(b+2)
                        
                    elif p%3 == 2 and (b+1) < max(features.location):
                        cod = ref[b-1]+ref[b]+ref[b+1]
                        m_cod = seq[b-1]+seq[b]+seq[b+1]                       
                        if ref[b-1] != seq[b-1]:
                            positions.append(b-1)
                        positions.append(b)
                        if ref[b+1] != seq[b+1]:
                            positions.append(b+1)

                    e = seq.id, positions, features.qualifiers['db_xref'][0], cod, translate(cod), m_cod, translate(m_cod), min(features.location), max(features.location)
                    result.append(e)
    if result:
        return result

def isIn(item, listIn):
    
    itemIn = item[0]
    if len(listIn) == 0:
        return False

    for it in listIn[0]:
        bOut = True

        # Seq name
        bOut = bOut and itemIn[0] == it[0]    

        # Position
        i = 0
        for pos in itemIn[1]:
            bOut = bOut and pos == itemIn[1][i]
            i += 1
        
        # GeneID
        bOut = bOut and itemIn[2] == it[2]

        # Ref
        bOut = bOut and itemIn[3] == it[3]

        # Seq
        bOut = bOut and itemIn[5] == it[5]

        # Min
        bOut = bOut and itemIn[7] == it[7]

        # Max
        bOut = bOut and itemIn[8] == it[8]

        if bOut:
            return True

    return False

def main():
    s = ""
    r = ""
    k = []
    l=0
    for f in glob.glob('Results/Kalign/*.fasta'):
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
                        diffs = []

                        # Output delle differenze, per ora solo in print
                        # TODO Output csv, txt, altro?

                        for position in range(j - l, j+1) :
                            if aln[i].seq[position].upper() != aln[0].seq[position].upper():
                                diff = find_genes(position, aln[0], aln[i])
                                if diff and not isIn(diff, diffs):
                                    diffs.append(find_genes(position, aln[0], aln[i]))
                        print(element)
                        if diffs:                     
                            print(diffs)
                            print("")
                        else:
                            print("Not found!")
                            print("")                            
                        k.append(element)
                        s = ""
                        r = ""
                        l = 0


if __name__ == "__main__":
	main()



