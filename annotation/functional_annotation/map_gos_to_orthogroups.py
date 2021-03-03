import os
import subprocess

def go_parents(target_go):
    if not os.path.exists("../go_hiers/hier_%s.txt" % target_go):
        cmd = ["python", "/Genomics/kocherlab/berubin/local/src/goatools/goatools/scripts/wr_hier.py", target_go, "-o", "../go_hiers/hier_%s.txt" % target_go, "--up", "--no_indent"]
        subprocess.call(cmd)
    reader = open("../go_hiers/hier_%s.txt" % target_go, 'rU')
    parents = []
    for line in reader:
        if line.startswith("*"):
            continue
        cur_line = line.strip()
        cur_line = cur_line.split()
        cur_parent = cur_line[0]
        parents.append(cur_parent)
    return parents


def add_parents(go_list, go_parents_dic):
    parent_go_list = []
    for go_num in go_list:
        parent_go_list.append(go_num)
    for go_num in go_list:
        if go_parents_dic.get(go_num, None) == None:
            go_parents_dic[go_num] = go_parents(go_num)
        go_parents_dic[go_num] = go_parents(go_num)
        for parent in go_parents_dic[go_num]:
            if parent not in parent_go_list:
                parent_go_list.append(parent)
    return parent_go_list, go_parents_dic
#    return parent_go_list, {}


og_go_counts = {}
og_fly_gene = {}
og_bee_gene = {}
og_odb_gene = {}
og_gene_counts = {}
go_parents_dic = {}
for species in ["AAUR", "APUR", "AVIR", "HLIG", "HQUA", "HRUB", "LALB", "LCAL", "LFIG", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "NMEL", "MGEN", "DNOV"]:
    print species
    line_count = 0
    reader = open("../%s/%s_trinotate_odb10.txt" % (species, species), 'rU')
    species_go_counts = {}
    species_fly_gene = {}
    species_bee_gene = {}
    species_odb_gene = {}
    for line in reader:
        line_count += 1
        if line_count % 1000 == 0:
#            break
            print line_count
        if line.startswith("#"):
            continue
        cur_line = line.split("\t")
        if cur_line[0] == "NA":
            continue
        cur_og = cur_line[0]
        cur_odb = cur_line[1]
        cur_fly = cur_line[2].split(";")[0]
        cur_bee = cur_line[3]
        if cur_line[4] == "NA":
            flygos = []
        else:
            flygos = cur_line[4].split(",")
        if cur_line[5] == "NA":
            beegos = []
        else:
            beegos = cur_line[5].split(",")
        go_list = []
        for field in cur_line[6:]:
            if field.startswith("GO:"):
                gos = field.split("`")
                for go_num in gos:
                    go_list.append(go_num.split("^")[0])
        go_list = list(set(go_list + flygos + beegos))
        gos, go_parents_dic = add_parents(go_list, go_parents_dic)
        go_list = list(set(gos))
        if og_go_counts.get(cur_og, None) == None:
            og_go_counts[cur_og] = []
            og_fly_gene[cur_og] = []
            og_bee_gene[cur_og] = []
            og_odb_gene[cur_og] = []
            og_gene_counts[cur_og] = 0.0
        if species_go_counts.get(cur_og, None) == None:
            species_go_counts[cur_og] = []
            species_fly_gene[cur_og] = []
            species_bee_gene[cur_og] = []
            species_odb_gene[cur_og] = []
#        og_go_counts[cur_og] = og_go_counts[cur_og] + go_list
        species_go_counts[cur_og] = species_go_counts[cur_og] + go_list
        species_fly_gene[cur_og].append(cur_fly)
        species_bee_gene[cur_og].append(cur_bee)
        species_odb_gene[cur_og].append(cur_odb)
#        og_fly_gene[cur_og].append(cur_fly)
#        og_bee_gene[cur_og].append(cur_bee)
#        og_odb_gene[cur_og].append(cur_odb)
#        og_gene_counts[cur_og] += 1
    for og_num, go_list in species_go_counts.items():
        og_go_counts[og_num] = og_go_counts[og_num] + list(set(go_list))
        og_fly_gene[og_num] = og_fly_gene[og_num] + list(set(species_fly_gene[og_num]))
        og_bee_gene[og_num] = og_bee_gene[og_num] + list(set(species_bee_gene[og_num]))
        og_odb_gene[og_num] = og_odb_gene[og_num] + list(set(species_odb_gene[og_num]))
        og_gene_counts[og_num] += 1

og_go_assigns = {}
og_fly_assigns = {}
og_bee_assigns = {}
og_odb_assigns = {}

for og_num, go_list in og_go_counts.items():
    for go_term in go_list:
        if go_list.count(go_term) / og_gene_counts[og_num] > 0.3:
            if og_go_assigns.get(og_num, None) == None:
                og_go_assigns[og_num] = []
            og_go_assigns[og_num].append(go_term)
    if og_go_assigns.get(og_num, None) != None:
        og_go_assigns[og_num] = list(set(og_go_assigns[og_num]))

outfile = open("og_go_assignments.txt", 'w')
for og_num, go_list in og_go_assigns.items():
    outfile.write("%s\t%s\n" % (og_num, ";".join(sorted(go_list))))
outfile.close()

for og_num, fly_list in og_fly_gene.items():
    og_fly_assigns[og_num] = "NA"
    for fly in fly_list:
        if fly_list.count(fly) / og_gene_counts[og_num] > 0.5:
            og_fly_assigns[og_num] = fly
            break

for og_num, bee_list in og_bee_gene.items():
    og_bee_assigns[og_num] = "NA"
    for bee in bee_list:
        if bee_list.count(bee) / og_gene_counts[og_num] > 0.5:
            og_bee_assigns[og_num] = bee
            break
        

for og_num, odb_list in og_odb_gene.items():
    og_odb_assigns[og_num] = "NA"
    for odb in odb_list:
        if odb_list.count(odb) / og_gene_counts[og_num] > 0.5:
            og_odb_assigns[og_num] = odb
            break
        
outfile = open("og_orthologs.txt", 'w')
for og_num, fly_gene in og_fly_assigns.items():
    outfile.write("%s\t%s\t%s\t%s\n" % (og_num, fly_gene, og_bee_assigns[og_num], og_odb_assigns[og_num]))
outfile.close()
