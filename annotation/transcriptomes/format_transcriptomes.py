from Bio import SeqIO
from glob import glob
import sys

species = sys.argv[1]

denovos = open("%s/%s_denovo.fasta" % (species, species), 'w')
ggs = open("%s/%s_gg.fasta" % (species, species), 'w')
denovo_counter = 0
gg_counter = 0
for transcripts in glob("%s/*Trinity*fasta" % species):
    print transcripts
    if not "-GG" in transcripts:
        reader = SeqIO.parse(transcripts, format = 'fasta')
        for rec in reader:
            denovos.write(">%s\n%s\n" % (rec.id.replace("TRINITY", "TRINITY" + str(denovo_counter)), str(rec.seq)))
        denovo_counter += 1
    else:
        reader = SeqIO.parse(transcripts, format = 'fasta')
        for rec in reader:
            ggs.write(">%s\n%s\n" % (rec.id.replace("TRINITY_GG", "TRINITY_GG" + str(gg_counter)), str(rec.seq)))
        gg_counter += 1
denovos.close()
ggs.close()
