import sys
from Bio import SeqIO

ingff = sys.argv[1]
longest_prots = sys.argv[2]
outgff = sys.argv[3]

long_ids = []
reader = SeqIO.parse(longest_prots, format = 'fasta')
for rec in reader:
    long_ids.append(rec.id)

outfile = open(outgff, 'w')
reader = open(ingff, 'rU')

for line in reader:
    if line.startswith("#"):
        continue
    cur_line = line.split()
    if cur_line[2] == "gene":
        outfile.write(line)
    elif cur_line[2] == "mRNA":
        cur_name = cur_line[8].split("Name=")[1][:-1]
        if cur_name in long_ids:
            outfile.write(line)
    else:
        parents = cur_line[8].split("Parent=")[1].split(";")[0]
        parent_list = parents.split(",")
        for cur_parent in parent_list:
            if cur_parent in long_ids:
                outfile.write(line.replace("Parent=%s;" % parents, "Parent=%s;" % cur_parent))
                break
outfile.close()
