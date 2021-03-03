from Bio import SeqIO
import sys

genome_file = sys.argv[1]
bacteria_blast = sys.argv[2]
hits_out = sys.argv[3]
if len(sys.argv) > 4:
    animal_blast = sys.argv[4]
    scafs_for_removal = sys.argv[5]
    remove_scafs = True
else:
    animal_blast = False
    remove_scafs = False

def overlap(mytuple, tuplelist):
    included = False
    i = 0
    while i < len(tuplelist):
        cur_tuple = tuplelist[i]
        if mytuple[0] <= cur_tuple[0] and mytuple[1] >= cur_tuple[0]-10:
            if mytuple[1] >= cur_tuple[1]:
                new_tuple = mytuple
                mytuple = new_tuple
                del tuplelist[i]
                i = -1
            else:
                new_tuple = (mytuple[0], cur_tuple[1])
                mytuple = new_tuple
                del tuplelist[i]
                i = -1
        elif mytuple[1] >= cur_tuple[1] and mytuple[0] <= cur_tuple[1] + 10:
            if mytuple[0] >= cur_tuple[0]:
                new_tuple = (cur_tuple[0], mytuple[1])
                mytuple = new_tuple
                del tuplelist[i]
                i = -1   
        elif mytuple[0] >= cur_tuple[0] and mytuple[1] <= cur_tuple[1]:
            included = True
            break
        i += 1
    if not included:
        tuplelist.append(mytuple)
    return tuplelist
        
def sum_tuples(tuplelist):
    tot_len = 0
    for t in tuplelist:
        tot_len += t[1] - t[0]
    return tot_len

alb_dic = {}

used_regions = {}
reader = SeqIO.parse(genome_file, format = 'fasta')
for rec in reader:
    alb_dic[rec.id] = rec.seq.tostring()
    used_regions[rec.id] = []

animal_hits = {}
if animal_blast:
    reader = open(animal_blast, 'rU')
    for line in reader:
        animal_hits[line.split()[0]] = float(line.split()[2])

reader = open(bacteria_blast, 'rU')
outfile = open(hits_out, 'w')
if remove_scafs:
    removefile = open(scafs_for_removal, 'w')

for line in reader:
    cur_line = line.split()
    if float(cur_line[3]) > 100:
        cur_scaff = cur_line[0]
        cur_start = int(cur_line[6])
        cur_end = int(cur_line[7])
        seq_name = "%s_%s_%s" % (cur_line[0], cur_line[6], cur_line[7])
        if seq_name in animal_hits.keys():
            if animal_hits[seq_name] > float(cur_line[2]):
                continue
        used_regions[cur_scaff] = overlap((cur_start, cur_end), used_regions[cur_scaff])

whole_scaffs = []
tot_len = 0
for k, v in used_regions.items():
    if len(v) > 0:
        if sum_tuples(v) > (len(alb_dic[k]) - alb_dic[k].count("N")) / 4.0 or (sum_tuples(v) > (len(alb_dic[k]) - alb_dic[k].count("N")) / 10.0 and len(v) > 10):
            whole_scaffs.append(k)
            if remove_scafs:
                removefile.write(k + "\n")
            outfile.write(">%s\n%s\n" % (k, alb_dic[k]))
            tot_len += len(alb_dic[k])
        else:
            print "%s\t%s\t%s" % (k, sorted(v), len(alb_dic[k]))
            for frag in sorted(v):
                outfile.write(">%s_%s_%s\n%s\n" % (k, frag[0], frag[1], alb_dic[k][frag[0]:frag[1]]))
outfile.close()
if remove_scafs:
    removefile.close()
