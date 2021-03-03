from Bio import SeqIO
import sys

rm_file = sys.argv[1]
outfile = sys.argv[2]
uniprot_blast = sys.argv[3]
dmel_blast = sys.argv[4]

blast_file = open(uniprot_blast, 'rU')

sim_seq_list = []

for line in blast_file:
    seqid = line.split()[0]
    pident = float(line.split()[2])
    length = float(line.split()[3])
    bitscore = float(line.split()[5])
    qlen = float(line.split()[6])
    slen = float(line.split()[7])
    if bitscore >= 100 and seqid not in sim_seq_list:
        sim_seq_list.append(seqid)
    elif length / qlen >= 0.50 and pident > 50 and seqid not in sim_seq_list:
        sim_seq_list.append(seqid)
    elif length / slen >= 0.50 and pident > 50 and seqid not in sim_seq_list:
        sim_seq_list.append(seqid)

blast_file.close()

blast_file = open(dmel_blast, 'rU')
for line in blast_file:
    seqid = line.split()[0]
    pident = float(line.split()[2])
    length = float(line.split()[3])
    bitscore = float(line.split()[5])
    qlen = float(line.split()[6])
    slen = float(line.split()[7])
    if bitscore >= 100 and seqid not in sim_seq_list:
        sim_seq_list.append(seqid)
    elif length / qlen >= 0.50 and pident > 50 and seqid not in sim_seq_list:
        sim_seq_list.append(seqid)
    elif length / slen >= 0.50 and pident > 50 and seqid not in sim_seq_list:
        sim_seq_list.append(seqid)

blast_file.close()

seq_file = SeqIO.parse(rm_file, format = 'fasta')
outfile = open(outfile, 'w')
for rec in seq_file:
    if rec.id not in sim_seq_list and len(rec.seq.tostring()) >= 80:
        outfile.write(">%s\n%s\n" % (rec.id, rec.seq.tostring()))

outfile.close()
