import sys

blast_file = sys.argv[1]
tax_output = sys.argv[2]

reader = open(blast_file, 'rU')

tax_dic = {}
tax_file = open("taxonomy.txt", 'rU') #taxonomy.txt is just a list of bacterial genome sequence GenBank IDs and the corresponding species names so that taxonomy can be assigned to the contaminating bacterial sequences.
for line in tax_file:
    cur_line = line.strip(">").split(" ")
    tax_dic[cur_line[0]] = cur_line[1] + "_" + cur_line[2]

hit_dic = {}

for line in reader:
    cur_seq = line.split()[0]
    cur_hit = line.split()[1]
    cur_e = float(line.split()[-2])
    if cur_seq in hit_dic.keys():
        newhit = False
        for hitit in hit_dic[cur_seq]:
            if cur_hit == hitit[0]:
                newhit = True
        if newhit:
            hit_dic[cur_seq].append((cur_hit, cur_e))
    else:
        hit_dic[cur_seq] = []
        hit_dic[cur_seq].append((cur_hit, cur_e))

precise_dic = {}
for k, v in hit_dic.items():
    smallest = v[0][1]
    smallname = v[0][0]
    for hit_tuple in v:
        if hit_tuple[1] < smallest:
            smallest = hit_tuple[1]
            smallname = hit_tuple[0]
    precise_dic[k] = smallname

outfile = open(tax_output, 'w')

for k, v in precise_dic.items():
    if "_" not in k:
        outfile.write("%s\t%s\n" % (k, tax_dic[v]))

for k, v in precise_dic.items():
    if "_" in k:
        outfile.write("%s\t%s\n" % (k, tax_dic[v]))

outfile.close()
