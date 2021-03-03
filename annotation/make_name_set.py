from Bio import SeqIO
import sys

prot_file = sys.argv[1]
gff_file = sys.argv[2]
in_id_list = sys.argv[3]
out_id_list = sys.argv[4]

reader = SeqIO.parse(prot_file, format = 'fasta')

seq_ids = []
for rec in reader:
    seq_ids.append(rec.id)

outfile = open(out_id_list, 'w')

gff_dic = {}
reader = open(gff_file, 'rU')
for line in reader:
    if "#FASTA" in line:
        break
    if line.startswith("#"):
        continue
    cur_line = line.split()
    if cur_line[2] == "mRNA":
        cur_attr = cur_line[8].split(";")
        cur_id = cur_attr[0].split("ID=")[1]
        cur_name = cur_attr[2].split("Name=")[1]    
        gff_dic[cur_id] = cur_name

reader = open(in_id_list, 'rU')
for line in reader:
    cur_line = line.split()
    if cur_line[0] in gff_dic.keys():
        outfile.write("%s\t%s\n" % (gff_dic[cur_line[0]], cur_line[1]))
    else:
        outfile.write(line)
outfile.close()
