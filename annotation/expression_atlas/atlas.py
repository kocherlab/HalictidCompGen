import subprocess
from Bio import SeqIO
import sys

BODY_PART = sys.argv[1] #options are "ant", "head", "abd", "duf"

SOCIAL_LIST = ["AAUR", "LMAR", "LZEP", "LMAL", "LCAL", "LPAU", "HLIG"]
SOLITARY_LIST = ["APUR", "LLEU", "LFIG", "LVIE", "LOEN", "HQUA"]
ANC_SOL_LIST = ["DNOV", "NMEL", "AVIR"]
POLY_LIST = ["MGEN", "LALB", "HRUB"]


def salmon(trans_dir, outdir):
    reader = open("barcode_map_combined.csv", 'rU')
    species_list = []
    data_dir = "/Genomics/kocherlab/berubin/data/rnaseq/all_demult_rnaseq"
    official_dic = {"calceatum" : "LCAL", "leucozonium" : "LLEU", "malachurum" : "LMAL", "pauxillum" : "LPAU", "pura": "APUR", "ligatus": "HLIG", "figueresi": "LFIG", "marginatum": "LMAR", "oenotherae": "LOEN", "vierecki": "LVIE", "zephyrum": "LZEP", "virescens": "AVIR", "aurata" : "AAUR", "albipes": "LALB"}
    for species in official_dic.values():
        cds_file = open("%s/%s_cds.fna" % (trans_dir, species), 'w')
        
        reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release_v2.1/%s/%s_OGS_v2.1_longest_isoform.cds.fasta" % (species, species), format = 'fasta')
        for rec in reader:
            cds_file.write(">%s\n%s\n" % (rec.id, str(rec.seq)))
        cds_file.close()
        cmd = ["/Genomics/kocherlab/berubin/local/src/Salmon-latest_linux_x86_64/bin/salmon", "index", "-t", "%s/%s_cds.fna" % (trans_dir, species), "-i", "%s/%s_trans_index" % (trans_dir, species), "--type", "quasi", "-k", "31"]
        subprocess.call(cmd)
    reader = open("barcode_map_combined.csv", 'rU')
    for line in reader:
        cur_line = line.split()
        if cur_line[1] != "callum":
            if cur_line[0].split("_")[-1] != BODY_PART:
                continue
        else:
            if cur_line[2] == BODY_PART:
                if BODY_PART == "ant" and cur_line[3] in ["pura", "marginatum"]:
                    continue
            else:
                continue
        cur_sample = cur_line[0]
        cur_species = official_dic[cur_line[3]]
        
        if "callum" in line:
            cmd = ["/Genomics/kocherlab/berubin/local/src/Salmon-latest_linux_x86_64/bin/salmon", "quant", "-i", "%s/%s_trans_index" % (trans_dir, cur_species), "-l", "IU", "-p", "6", "--gcBias", "-1", "%s/%s_R1.fastq.gz" % (data_dir, cur_sample), "-2", "%s/%s_R2.fastq.gz" % (data_dir, cur_sample), "-o", "%s/%s_%s_salmon" % (outdir, cur_species, cur_sample)]
            subprocess.call(cmd)
        else:
            cmd = ["/Genomics/kocherlab/berubin/local/src/Salmon-latest_linux_x86_64/bin/salmon", "quant", "-i", "%s/%s_trans_index" % (trans_dir, cur_species), "-l", "ISR", "-p", "6", "--gcBias", "-1", "%s/%s_R1.fastq.gz" % (data_dir, cur_sample), "-2", "%s/%s_R2.fastq.gz" % (data_dir, cur_sample), "-o", "%s/%s_%s_salmon" % (outdir, cur_species, cur_sample)]
            subprocess.call(cmd)
        

def read_orthofinder():
    reader = open("/Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/protein_seqs/OrthoFinder/Results_Jan29_4/Orthogroups/Orthogroups.txt", 'rU')
    og_dic = {}
    for line in reader:
        cur_line = line.split()
        cur_og = cur_line[0][4:-1]
        cur_og = int(cur_og) + 10000
        if len(cur_line) == 2:
            continue
        for gene in cur_line[1:]:
            og_dic[gene] = "OG_" + str(cur_og)
    return og_dic


def salmon_output(outdir):
    og_dic = read_orthofinder()
    reader = open("barcode_map_combined.csv", 'rU')
    species_list = []
    official_dic = {"calceatum" : "LCAL", "leucozonium" : "LLEU", "malachurum" : "LMAL", "pauxillum" : "LPAU", "pura": "APUR", "ligatus": "HLIG", "figueresi": "LFIG", "marginatum": "LMAR", "oenotherae": "LOEN", "vierecki": "LVIE", "zephyrum": "LZEP", "virescens": "AVIR", "aurata" : "AAUR", "albipes": "LALB"}
    outfile = open("salmon_%s.txt" % BODY_PART, 'w')
    outfile.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\tSocial\tSpecies\n")
    for line in reader:
        cur_line = line.split()
        if cur_line[1] != "callum":
            if cur_line[0].split("_")[-1] != BODY_PART:
                continue
        else:
            if cur_line[2] == BODY_PART:
                if BODY_PART == "ant" and cur_line[3] in ["pura", "marginatum"]:
                    continue
            else:
                continue
        cur_sample = cur_line[0]
        cur_species = official_dic[cur_line[3]]
        if cur_species in species_list:
            continue
        species_list.append(cur_species)
        salmon = open("%s/%s_%s_salmon/quant.sf" % (outdir, cur_species, cur_sample), 'rU')

        if cur_species in SOCIAL_LIST:
            social_state = "soc"
        elif cur_species in SOLITARY_LIST:
            social_state = "sol"
        elif cur_species in POLY_LIST:
            social_state = "poly"
        else:
            social_state = "ancsol"
        chemo_dic = {}
        for gene_line in salmon:
            if gene_line.startswith("Name"):
                continue
            chemo_dic[gene_line.split()[0]] = gene_line.split()
        for gene, gene_line in chemo_dic.items():
            genlen = int(gene_line[1])
            numreads = float(gene_line[4])
            if gene in og_dic.keys():
                cur_og = og_dic[gene]
                del og_dic[gene]
            else:
                cur_og = "NA"
            outfile.write("\t".join(gene_line) + "\t" + "\t" + social_state + "\t" + cur_species + "\t" + cur_og + "\n")
    outfile.close()


def main():
    salmon("/Genomics/kocherlab/berubin/transcriptomes/expression/hic/good_transcripts", "/Genomics/kocherlab/berubin/transcriptomes/expression/hic/salmon_out")
    salmon_output("/Genomics/kocherlab/berubin/transcriptomes/expression/hic/salmon_out")
if __name__ == '__main__':
    main()
