import random
from glob import glob
from Bio import SeqIO
import copy
from ete3 import PhyloTree

def main():
    assign_genes()
    
def assign_genes():
    og_dic = {}
    reader = open("/Genomics/kocherlab/berubin/annotation/hic/trinotate/NMEL/NMEL_trinotate_odb10.txt", 'rU') #This file is just used as a map between the NMEL genes and OG numbers.
    for line in reader:
        cur_line = line.split()
        og_dic[cur_line[6][0:-3]] = cur_line[0]
    gff_dic = read_gff("/Genomics/kocherlab/berubin/comparative/halictids/genome_alignments/gffs/NMEL.gff") #This is just the OGS gff file
    outfile = open("cne_proximity.txt", 'w')
    outfile.write("NCAR\tOG\tNMEL_gene\tproxtype\tNMEL_locus\n")
    reader = open("/Genomics/kocherlab/berubin/comparative/halictids/halictid_selection/halictid_ncar_ortho.index", 'rU') #This is the index file produced by the first step of the genomics pipline. It is "halictid_cnee_ortho.index" on the google drive.
    ncar_list = []
    nmel_missing = []
    for line in reader:
        nmel_found = False
        cur_line = line.split()
        seq_list = cur_line[3].split(",")
        for seq in seq_list:
            if seq.startswith("NMEL"):
                ncar_list.append(seq)
                nmel_found = True
        if not nmel_found:
            nmel_missing.append(cur_line[1])
    print len(nmel_missing)
    coord_dic = {}
    proxim_dic = {}
    for ncar in ncar_list:
        cur_deets = ncar.split(":")
        cur_ce = cur_deets[5]
        cur_scaf = cur_deets[1]
        cur_start = int(cur_deets[2])
        cur_end = int(cur_deets[3])
        cur_middle = (cur_start + cur_end) / 2
        if cur_scaf not in coord_dic.keys():
            coord_dic[cur_scaf] = []
        cur_proxim = process_coord(cur_middle, cur_scaf, gff_dic)
        print ncar
        print cur_proxim
        proxim_dic[ncar] = cur_proxim
#        coord_dic[cur_scaf].append(cur_start)
        cur_species = []
#        for species in phylo.get_leaves():
#            cur_species.append(species.name)
        if proxim_dic[ncar][1] == "NA":
            outfile.write("%s\t%s\t%s\t%s\t%s\n" % (cur_ce, "NA", proxim_dic[ncar][1], proxim_dic[ncar][0], ncar))
        else:
            outfile.write("%s\t%s\t%s\t%s\t%s\n" % (cur_ce, og_dic[proxim_dic[ncar][1]], proxim_dic[ncar][1], proxim_dic[ncar][0], ncar))
        outfile.flush()
    outfile.close()
            
def process_coord(start_coord, scaf, gene_dic):
    if scaf not in gene_dic.keys():
        return "intergenic", "NA"
    added = False
    coord_proxim_dic = {}
    for gene in gene_dic[scaf]:
        cur_proxim = gene.proximity(start_coord)
        if not cur_proxim:
            continue
        coord_proxim_dic[gene.name] = cur_proxim
        added = True
    if added:
        best_proxim = choose_best_proxim(coord_proxim_dic)
        for gene in gene_dic[scaf]:
            if gene.name == best_proxim[0]:
                best_gene = gene
                break
        proxtype = best_proxim[1]
        return proxtype, best_gene.name
        added = True
#        break
    if not added:
        return "intergenic", "NA"

def process_coord_manyhits(start_coord, scaf, gene_dic):
    mod_name_dic = get_old_og_dic()
    if scaf not in gene_dic.keys():
        return [("intergenic", "NA", "NA", "NA", "NA")]
    added = False
    hits_list = []
    for gene in gene_dic[scaf]:
        cur_proxim = gene.proximity(start_coord)
        if not cur_proxim:
            continue
        express_bias = gene.expression_bias
        wq_express_bias = gene.wq_expression_bias
        proxtype = cur_proxim[0]
        if gene.name in mod_name_dic.keys():
            hits_list.append((proxtype, gene.name, mod_name_dic[gene.name], express_bias, wq_express_bias))
        else:
            hits_list.append((proxtype, gene.name, "NA", express_bias, wq_express_bias))
        added = True
#        break
    if not added:
        return [("intergenic", "NA", "NA", "NA", "NA")]
    else:
        return hits_list

def read_gff(gff_file):
    reader = open(gff_file, 'rU')
    ingene = False
    first_gene = True
    scaf_dic = {}
    for line in reader:
        cur_line = line.split()
        cur_scaf = cur_line[0]
        if cur_scaf not in scaf_dic.keys():
            scaf_dic[cur_scaf] = []
        if cur_line[2] == "gene":
            if not first_gene:
                cur_gene.get_introns()
                scaf_dic[cur_gene.scaf].append(cur_gene)
            first_gene = False
            cur_id = cur_line[-1].split("ID=")[1].split(";Name")[0]
            cur_gene = Gene(cur_id, int(cur_line[3]), int(cur_line[4]), cur_line[6], cur_line[0])
        elif cur_line[2] == "CDS":
            cur_gene.add_cds(int(cur_line[3]), int(cur_line[4]))
        elif cur_line[2] == "five_prime_UTR":
            cur_gene.add_five(int(cur_line[3]), int(cur_line[4]))
        elif cur_line[2] == "three_prime_UTR":
            cur_gene.add_three(int(cur_line[3]), int(cur_line[4]))
    scaf_dic[cur_scaf].append(cur_gene)
    return scaf_dic


def choose_best_proxim(coord_proxim_dic):
    
    proxtype_counts = {}
    proxtype_genes = {}
    for proxtype in ["intron", "three_utr", "five_utr", "promoter", "stream"]:
        
        proxtype_counts[proxtype] = []
    best_per_gene = {}
    for gene_name, proxim in coord_proxim_dic.items():
        if proxim[0] in ["downstream", "upstream"]:
            proxtype_counts["stream"].append((gene_name, proxim[0], proxim[1])) 
        else:
            proxtype_counts[proxim[0]].append((gene_name, proxim[0], proxim[1]))
    for proxtype in ["intron", "three_utr", "five_utr", "promoter", "stream"]:
#    for proxtype, prox_list in proxtype_counts.items():
        prox_list = proxtype_counts[proxtype]
        if len(prox_list) == 0:
            continue
        if proxtype == "intron":
            return random.choice(prox_list)
        elif proxtype == "three_utr":
            return random.choice(prox_list)
        elif proxtype == "five_utr":
            return random.choice(prox_list)
        elif proxtype == "promoter":
            return choose_closest(prox_list)
        else:
            return choose_closest(prox_list)

def choose_closest(prox_list):
    closest_coord = prox_list[0][2]
    closest_prox = prox_list[0]
    for prox in prox_list:
        if prox[2] < closest_coord:
            closest_coord = prox[2]
            closest_prox = prox
    return closest_prox


class Gene:
    def __init__(self, gene_name, gene_start, gene_end, gene_strand, gene_scaf):
        self.name = gene_name
        self.start = gene_start
        self.end = gene_end
        self.strand = gene_strand
        self.scaf = gene_scaf
        self.cds = []
        self.five_utrs = []
        self.three_utrs = []
        self.introns = []
        
    def add_cds(self, cur_start, cur_end):
        self.cds.append((cur_start+1, cur_end-1))

    def add_five(self, cur_start, cur_end):
        self.five_utrs.append((cur_start-1, cur_end+1))

    def add_three(self, cur_start, cur_end):
        self.three_utrs.append((cur_start-1, cur_end+1))

    def get_introns(self):
        intron_list = []
        for x in range(len(self.cds)-1):
            if self.strand == "+":
                cur_intron = (self.cds[x][1], self.cds[x+1][0])
            else:
                cur_intron = (self.cds[x+1][0], self.cds[x][1])
            intron_list.append(cur_intron)
        self.introns = intron_list
        

    def coord_in(self, target_coord, target_list):
        for coords in target_list:
#            if self.strand == "+":
            if target_coord >= coords[0] and target_coord <= coords[1]:
                return True
            else:
                if target_coord >= coords[1] and target_coord <= coords[0]:
                    return True

    def proximity(self, target_coord):
        if target_coord < self.end + 5000 and target_coord > self.start - 5000:

            if self.coord_in(target_coord, self.introns):
                return ("intron", -9)
            elif self.coord_in(target_coord, self.five_utrs):
                return ("five_utr", -9)
            elif self.coord_in(target_coord, self.three_utrs):
                return ("three_utr", -9)
            elif self.strand == "+":
                if target_coord <= self.start:
                    if target_coord >= self.start - 1500:
                        return ("promoter", self.start - target_coord)
                    else:
                        return ("upstream", self.start - target_coord)
                else:
                    return ("downstream", target_coord - self.end)
            elif self.strand == "-":
                if target_coord >= self.end:
                    if target_coord <= self.end + 1500:
                        return ("promoter", target_coord - self.end)
                    else:
                        return ("upstream", target_coord - self.end)
                else:
                    return ("downstream", self.start - target_coord)

            elif self.coord_in(target_coord, self.cds):
                print self.name
                print target_coord
                return ("cds", -9)
        else:
            return False

if __name__ == '__main__':
    main()



