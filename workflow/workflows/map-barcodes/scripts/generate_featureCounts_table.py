# snakemake workflow version

# input[0]: a sam file with barcodes and gene identifiers (output of featureCounts)
# output[0]: a table

from common.utils import chunks

# Expecting these columns:
# 1 	QNAME 	String 	Query template NAME
# 2 	FLAG 	Int 	bitwise FLAG
# 3 	RNAME 	String 	References sequence NAME
# 4 	POS 	Int 	1- based leftmost mapping POSition
# 5 	MAPQ 	Int 	MAPping Quality
# 6 	CIGAR 	String 	CIGAR String
# 7 	RNEXT 	String 	Ref. name of the mate/next read
# 8 	PNEXT 	Int 	Position of the mate/next read
# 9 	TLEN 	Int 	observed Template LENgth
# 10 	SEQ 	String 	segment SEQuence
# 11 	QUAL 	String 	ASCII of Phred-scaled base QUALity+33
# 12    NH
# 13    HI
# 14    AS
# 15    nM
# 16    bc      String  Cell barcode
# 17    um      String  UMI
# 18    XS
# 19    XN
# 20    XT      String  Gene identifier

print("Writing output file: {}".format(snakemake.output[0]))
fi = open(snakemake.input[0], 'r') 
fo = open(snakemake.output[0], 'w', buffering=32768)
fo.write('\t'.join(['readid','reference','position','mapquality','cellbc','umi','xs','xn','geneid']) + '\n')
for line in fi:
    if not line.startswith('@'):
        line = line[:-1] # removing '\n' at end of line
        e = line.split('\t')
        ne = len(e)-1
        bc = e[15].split(':')[2]
        um = e[16].split(':')[2]
        xs = e[17].split(':')[2]
        if (ne>17):
            xn = e[18].split(':')[2]
            xt = e[19].split(':')[2]
        else:
            xn = 'NA'
            xt = 'NA'
        fo.write('\t'.join([e[0],e[2],e[3],e[4],bc,um,xs,xn,xt]) + '\n')
fi.close()
fo.close() 
