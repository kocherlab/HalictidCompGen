import numpy

def check_paralog(gene_list):
    species_present = []
    for gene in gene_list:
        cur_species = gene[0:4]
        if cur_species in species_present:
            return True
        species_present.append(cur_species)
    return False

def get_paralog_species(gene_list):
    species_present = []
    paralog_species = []
    for gene in gene_list:
        cur_species = gene[0:4]
        if cur_species in species_present:
            paralog_species.append(cur_species)
        species_present.append(cur_species)
    return paralog_species


def count_species(gene_list):
    species_present = []
    for gene in gene_list:
#        cur_species = gene.split("_")[0]
        cur_species = gene[0:4]
        species_present.append(cur_species)
    return len(list(set(species_present)))

def check_identical(oglist1, oglist2):
    for og1 in oglist1:
        if og1 not in oglist2:
            return False
    for og2 in oglist2:
        if og2 not in oglist1:
            return False
    return True

def check_subset_of_2(oglist1, oglist2):
    for og1 in oglist1:
        if og1 not in oglist2:
            return False
    return True

def same_species(oglist1, oglist2):
    speclist1 = []
    for og1 in oglist1:
        speclist1.append(og1[0:4])
    speclist2 = []
    for og2 in oglist2:
        speclist2.append(og2[0:4])
    speclist1 = list(set(speclist1))
    speclist2 = list(set(speclist2))
    return check_identical(speclist1, speclist2)

#reader = open("/Genomics/kocherlab/berubin/comparative/halictids/orthology/oma/OMA_Output_v1/OrthologousGroups_renamed.txt", 'rU')
reader = open("/Genomics/kocherlab/berubin/comparative/halictids/orthology/oma/OMA_Output_v1/HOG_groups.txt", 'rU')
oma_dic = {}
oma_species_dic = {}
gene_oma_dic = {}
oma_by_species = {}
for line in reader:
    if line.startswith("#"):
        continue
    cur_line = line.strip().split()
    cur_og = cur_line[0]
    gene_list = []
    species_list = []
    for gene in cur_line[1:]:
        cur_gene = gene.split(":")[1]
        cur_species = cur_gene[0:4]
        gene_list.append(cur_gene)
        species_list.append(cur_species)
        gene_oma_dic[cur_gene] = cur_og
    oma_dic[cur_og] = gene_list
    oma_species_dic[cur_og] = species_list
    for gene in gene_list:
        species = gene[0:4]
        if species not in oma_by_species.keys():
            oma_by_species[species] = {}
        oma_by_species[species][gene] = gene_list

reader = open("/Genomics/kocherlab/berubin/annotation/hic/orthomcl/groups.txt", 'rU')

fullmcl_dic = {}
fullmcl_species_dic = {}
fullmcl_by_species = {}
for line in reader:
    cur_line = line.split()
    gene_list = []
    species_list = []
    for gene in cur_line[1:]:
        gene_list.append(gene.split("|")[1])
        species_list.append(gene.split("|")[0])
    fullmcl_dic[cur_line[0][0:-1]] = gene_list
    fullmcl_species_dic[cur_line[0][0:-1]] = list(set(species_list))
    for gene in gene_list:
        species = gene[0:4]
        if species not in fullmcl_by_species.keys():
            fullmcl_by_species[species] = {}
        fullmcl_by_species[species][gene] = gene_list

reader = open("/Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/protein_seqs/OrthoFinder/Results_Jan29_4/Orthogroups/Orthogroups.txt", 'rU')

orthof_dic = {}
orthof_species_dic = {}
orthof_by_species = {}
for line in reader:
    cur_line = line.split()
    gene_list = []
    species_list = []
    for gene in cur_line[1:]:
        gene_list.append(gene)
        species_list.append(gene[0:4])
    orthof_dic[cur_line[0][0:-1]] = gene_list
    orthof_species_dic[cur_line[0][0:-1]] = list(set(species_list))
    for gene in gene_list:
        species = gene[0:4]
        if species not in orthof_by_species.keys():
            orthof_by_species[species] = {}
        orthof_by_species[species][gene] = gene_list

mcl_dic = {}
fly_dic = {}
bee_dic = {}
bee_species_dic = {}
species_genes = {}
for species in ["AAUR", "APUR", "AVIR", "HLIG", "HQUA", "HRUB", "LALB", "LCAL", "LFIG", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "NMEL"]:
    bee_species_dic[species] = {}
    species_genes[species] = []
    mcl_count = 0
    fly_count = 0
    bee_count = 0
    mcl_fly_count = 0
    mcl_bee_count = 0
    fly_bee_count = 0
    mcl_fly_bee_count = 0
    gene_count = 0
    mcl = False
    bee = False
    fly = False
    reader = open("%s/%s_trinotate_odb10.txt" % (species, species), 'rU')
    for line in reader:
        if line.startswith("#OG"):
            continue
        gene_count += 1
        cur_line = line.split()
        cur_gene = cur_line[7]
        species_genes[species].append(cur_gene)
        if cur_line[1] != "NA":
            bee = True
            bee_count += 1
            if bee_dic.get(cur_line[1], None) == None:
                bee_dic[cur_line[1]] = []
            bee_dic[cur_line[1]].append(cur_gene)
            bee_species_dic[species][cur_gene] = cur_line[1]

def get_method_counts(og_dic, species_genes, odb_dic):
    mcl_count = 0
    odb_count = 0
    mcl4_odb_count = 0
    mcl10_odb_count = 0
    mcl10_count = 0
    mcl_species10 = 0
    mcl_species10_thisparalog = 0
    mcl_species10_anyparalog = 0
    for gene in species_genes:
        if gene in odb_dic.keys():
            odb_count += 1
        gene_list = og_dic.get(gene, None)
        if gene_list == None:
            continue
#        for mcl, gene_list in og_dic.items():
        if gene in gene_list:
            mcl_count += 1
            if count_species(gene_list) >= 4:
                if gene in odb_dic.keys():
                    mcl4_odb_count += 1
            if len(gene_list) >= 10:
                mcl10_count += 1
                if count_species(gene_list) >= 10:
                    if gene in odb_dic.keys():
                        mcl10_odb_count += 1
                    mcl_species10 += 1
                    if gene[0:4] in get_paralog_species(gene_list):
                        mcl_species10_thisparalog += 1
                    if len(get_paralog_species(gene_list)) > 0:
                        mcl_species10_anyparalog += 1

    return [mcl_count, odb_count, mcl10_count, mcl_species10, mcl_species10_thisparalog, mcl_species10_anyparalog, mcl4_odb_count, mcl10_odb_count]

def write_species_summary():
    outfile = open("assess_methods_counts.txt", 'w')
    outfile.write("species\t#genes\t#mcls\t#odb\t#mcl10\tmclsp10\tmclsp10thispara\tmclsp10anypara\tmcl4odb\tmcl10odb\t#orfs\t#odb\t#orf10\torfsp10\torfsp10thispara\torfsp10anypara\torf4odb\torf10odb\t#omas\t#odb\t#oma10\tomasp10\tomasp10thispara\tomasp10anypara\toma4odb\toma10odb\n")
    
    for species in ["AAUR", "APUR", "AVIR", "HLIG", "HQUA", "HRUB", "LALB", "LCAL", "LFIG", "LLEU", "LMAL", "LMAR", "LOEN", "LPAU", "LVIE", "LZEP", "NMEL"]:
        print species
        mcl_deets = get_method_counts(fullmcl_by_species[species], species_genes[species], bee_species_dic[species])
        orthof_deets = get_method_counts(orthof_by_species[species], species_genes[species], bee_species_dic[species])
        oma_deets = get_method_counts(oma_by_species[species], species_genes[species], bee_species_dic[species])
        outline = list(map(str, mcl_deets + orthof_deets + oma_deets))
        outfile.write("%s\t%s\t%s\n" % (species, len(species_genes[species]), "\t".join(outline)))
        outfile.flush()
    outfile.close()

#write_species_summary()

def get_stats(ortho_dic, min_size):
    mcl10 = 0
    mcl_paralogs = 0
    mcl_sizes = []
    for mcl, gene_list in ortho_dic.items():
        if count_species(gene_list) >= min_size:
            mcl10 += 1
            mcl_sizes.append(count_species(gene_list))
            if check_paralog(gene_list):
                mcl_paralogs += 1
    percents = "\t".join(list(map(str, numpy.percentile(mcl_sizes, [10, 20, 30, 40, 50]).tolist())))
    return numpy.median(mcl_sizes), round(numpy.mean(mcl_sizes), 3), numpy.sum(mcl_sizes), percents, mcl10, mcl_paralogs          
 
outfile = open("method_stats.txt", 'w')
outfile.write("Method\tminspec\tmedsiz\tmeansiz\ttotgen\tperc10\tperc20\tperc30\tperc40\tperc50\tnumogs\togsw/para\n")
for minsize in [4, 10, 15]:
    print minsize
    outline = list(map(str, get_stats(fullmcl_dic, minsize)))
    outfile.write("MCL\t%s\t%s\n" % (minsize, "\t".join(outline)))
    outline = list(map(str, get_stats(orthof_dic, minsize)))
    outfile.write("OrthoF\t%s\t%s\n" % (minsize, "\t".join(outline)))
    outline = list(map(str, get_stats(oma_dic, minsize)))
    outfile.write("OMA\t%s\t%s\n" % (minsize, "\t".join(outline)))
outfile.close()


def compare_og_dics(ogdic1, ogdic2):
    match_count = 0
    for mclogs in ogdic1.values():
        if count_species(mclogs) < 10:
            continue
        matcher = False
        sub_match = False
        for omaogs in ogdic2.values():
            if check_identical(mclogs, omaogs):
                matcher = True
                break
        if matcher:
            match_count += 1
    return match_count

def compare_dic_odb(target_dic, odb_dic):
    new_dic = {}
    for ognum, gene_list in target_dic.items():
        new_genes = []
        for gene in gene_list:
            if gene[0:4] in ["Dnov", "MGEN", "DNOV", "Mgen"]:
                continue
            new_genes.append(gene)
        new_dic[ognum] = new_genes
    return compare_og_dics(new_dic, odb_dic)

print "Exacts between mcl and oma: %s" % compare_og_dics(fullmcl_dic, oma_dic)
print "Exacts between mcl and orthofinder: %s" % compare_og_dics(fullmcl_dic, orthof_dic)
print "Exacts between oma and orthofinder: %s" % compare_og_dics(oma_dic, orthof_dic)
print "Exacts between mcl and odb: %s" % compare_dic_odb(fullmcl_dic, bee_dic)
print "Exacts between oma and odb: %s" %compare_dic_odb(oma_dic, bee_dic)

print "Exacts between orthofinder and odb: %s" %compare_dic_odb(orthof_dic, bee_dic)
