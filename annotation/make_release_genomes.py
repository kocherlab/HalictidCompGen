from Bio import SeqIO
import sys
import os

inspecies = sys.argv[1]

if not os.path.exists("/Genomics/kocherlab/berubin/official_release_v2.1/%s" % inspecies):
    os.mkdir("/Genomics/kocherlab/berubin/official_release_v2.1/%s" % inspecies)

reader = SeqIO.parse("/Genomics/kocherlab/berubin/assembly/hic/%s/%s_genome_v2.0.fasta" % (inspecies, inspecies), format = 'fasta')

seq_dic = {}
for rec in reader:
    seq_dic[rec.id] = str(rec.seq)

reader = open("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v2.1.gff3" % (inspecies, inspecies), 'rU')

gff_dic = {}
for line in reader:
    cur_line = line.strip().split("\t", 1)
    cur_scaf = cur_line[0]
    if cur_scaf not in gff_dic.keys():
        gff_dic[cur_scaf] = []
    gff_dic[cur_scaf].append(cur_line[1])

rename_dic = {}
if os.path.exists("/Genomics/kocherlab/berubin/assembly/hic/%s/%s_renames.txt" % (inspecies, inspecies)):
    reader = open("/Genomics/kocherlab/berubin/assembly/hic/%s/%s_renames.txt" % (inspecies, inspecies), 'rU')
    for line in reader:
        if line.startswith("OldName"):
            continue
        cur_line = line.split()
        seq_dic[cur_line[1]] = seq_dic[cur_line[0]]
        del seq_dic[cur_line[0]]
        gff_dic[cur_line[1]] = gff_dic[cur_line[0]]
        del gff_dic[cur_line[0]]
        rename_dic[cur_line[0]] = cur_line[1]

removal_genes = []
if os.path.exists("/Genomics/kocherlab/berubin/assembly/hic/%s/%s_other_taxon.txt" % (inspecies, inspecies)):
    discard_list = []
    reader = open("/Genomics/kocherlab/berubin/assembly/hic/%s/%s_other_taxon.txt" % (inspecies, inspecies), 'rU')
    for line in reader:
        discard_list.append(int(line.strip()))
    for k, v in seq_dic.items():
        cur_index = int(k.split("_")[-1])
        if cur_index in discard_list:
            del seq_dic[k]
            if k in gff_dic.keys():
                for gff_line in gff_dic[k]:
                    if "gene" in gff_line:
                        gene_name = gff_line.split(";")[1].split("Name=")[1]
                        removal_genes.append(gene_name)
                del gff_dic[k]

if os.path.exists("/Genomics/kocherlab/berubin/assembly/hic/%s/%s_alternatives.txt" % (inspecies, inspecies)):
    reader = open("/Genomics/kocherlab/berubin/assembly/hic/%s/%s_alternatives.txt" % (inspecies, inspecies), 'rU')
    for line in reader:
        if line.startswith("Alt"):
            continue
        cur_line = line.split()
        print cur_line
        if "%s_scaf_%s" % (inspecies, cur_line[0]) in seq_dic.keys():
            cur_scaf = "%s_scaf_%s" % (inspecies, cur_line[0])
        elif "%s_chr_%s" % (inspecies, cur_line[0]) in seq_dic.keys():
            cur_scaf = "%s_chr_%s" % (inspecies, cur_line[0])
        new_scaf = cur_scaf + "_altto_%s" % cur_line[1]
        seq_dic[new_scaf] = seq_dic[cur_scaf]
        del seq_dic[cur_scaf]
        if cur_scaf in gff_dic.keys():
            for gff_line in gff_dic[cur_scaf]:
                if "gene" in gff_line:
                    gene_name = gff_line.split(";")[1].split("Name=")[1]
                    removal_genes.append(gene_name)
            del gff_dic[cur_scaf]

gff_good_nums = []
outgff = open("/Genomics/kocherlab/berubin/official_release_v2.1/%s/%s_OGS_v2.1.gff3" % (inspecies, inspecies), 'w')
for scaf, gff_lines in gff_dic.items():
    new_scaf = scaf.replace("scaf", "unplaced")
    gff_good_nums.append(int(new_scaf.split("_")[-1]))
    for line in gff_lines:
        outgff.write(new_scaf + "\t" + line + "\n")
outgff.close()

outgff = open("/Genomics/kocherlab/berubin/official_release_v2.1/%s/%s_OGS_v2.1_longest_isoform.gff3" % (inspecies, inspecies), 'w')
ingff = open("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v2.1_longest_isoform.gff3" % (inspecies, inspecies), 'rU')
for line in ingff:
    cur_line = line.split()
    cur_scaf = cur_line[0]
    if cur_scaf in gff_dic.keys():
        outgff.write(line.replace("scaf", "unplaced"))
    elif cur_scaf in rename_dic.keys():
        if int(cur_scaf.split("_")[-1]) in gff_good_nums:
            outgff.write(line.replace(cur_scaf, rename_dic[cur_scaf]).replace("scaf", "unplaced"))
outgff.close()

                
outgenome = open("/Genomics/kocherlab/berubin/official_release_v2.1/%s/%s_genome_v2.1.fasta" % (inspecies, inspecies), 'w')
for scaf, seq in seq_dic.items():
    new_scaf = scaf.replace("scaf", "unplaced")
    outgenome.write(">%s\n" % new_scaf)
    prev_index = 0
    cur_index = 80
    while prev_index < len(seq):
        outgenome.write("%s\n" % (seq[prev_index:cur_index]))
        prev_index = cur_index
        cur_index = cur_index + 80
outgenome.close()


for seq_file in ["longest_isoform.cds", "longest_isoform.gene", "longest_isoform.pep", "longest_isoform.trans", "longest_isoform.upstream3000", "pep", "trans"]:
    outfile = open("/Genomics/kocherlab/berubin/official_release_v2.1/%s/%s_OGS_v2.1_%s.fasta" % (inspecies, inspecies, seq_file), 'w')
    reader = SeqIO.parse("/Genomics/kocherlab/berubin/official_release/%s/%s_OGS_v2.1_%s.fasta" % (inspecies, inspecies, seq_file), format = 'fasta')
    for rec in reader:

        if rec.id[0:10] not in removal_genes:
            outfile.write(">%s\n%s\n" % (rec.id, str(rec.seq)))
outfile.close()
        
