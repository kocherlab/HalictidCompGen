from Bio import SeqIO
import sys, os

ingenome = sys.argv[1]
inprots = sys.argv[2]
intranscripts = sys.argv[3]
ingff = sys.argv[4]
bac_scafs = sys.argv[5]
outdir = sys.argv[6]
outprefix = sys.argv[7]

if not os.path.isdir("%s/%s" % (outdir, outprefix)):
    os.mkdir("%s/%s" % (outdir, outprefix))

bac_list = []
reader = open(bac_scafs, 'rU')
for line in reader:
    bac_list.append(line.split()[0])

outgenome = open("%s/%s/%s_genome_v2.1.fasta" % (outdir, outprefix, outprefix), 'w')
reader = SeqIO.parse(ingenome, format = 'fasta')
scaf_count = 0
scaf_renames = {}
for rec in reader:
    if rec.id not in bac_list:
        outgenome.write(">%s\n" % (rec.id))
        prev_index = 0
        cur_index = 80
        while prev_index < len(rec.seq):
            outgenome.write("%s\n" % (str(rec.seq)[prev_index:cur_index]))
            prev_index = cur_index
            cur_index = cur_index + 80
#        outgenome.write(">%s_scaf_%s\n%s\n" % (outprefix, scaf_count, str(rec.seq)))
        scaf_renames[rec.id] = rec.id 
        scaf_count += 1
outgenome.close()

bac_genes = []
outgff = open("%s/%s/%s_OGS_v2.1.gff3" % (outdir, outprefix, outprefix), 'w')
reader = open(ingff, 'rU')
for line in reader:
    if "#FASTA" in line:
        break
    if line.startswith("#"):
        continue
    cur_line = line.split()

    if cur_line[0] in bac_list:
        if cur_line[2] == "gene":
            cur_id = cur_line[8].split(";")[0].split("ID=")[1]
            bac_genes.append(cur_id)
        continue
    new_scaf = scaf_renames[cur_line[0]]

        
    if "Alias" in cur_line[8]:
        if cur_line[2] == "mRNA":
            cur_attr = cur_line[8].split(";")
            cur_id = cur_attr[0].split("ID=")[1]
            cur_name = cur_attr[2].split("Name=")[1]
            line = line.replace(cur_name, cur_id)
        outgff.write(line.split("Alias")[0].replace(cur_line[0], new_scaf) + "\n")
    else:
        outgff.write(line.replace(cur_line[0], new_scaf))
outgff.close()

def longest_iso(gene_list, iso_dic):
    longest_list = []
    for gene in gene_list:
        cur_isos = []
        for k, v in iso_dic.items():
            cur_id = k.split("-")[0]
            if cur_id == gene:
                cur_isos.append(k)
        longest = ""
        longest_length = 0
        for iso in cur_isos:
            if len(iso_dic[iso]) > longest_length:
                longest = iso
                longest_length = len(iso_dic[iso])
        
        longest_list.append(longest)
    return longest_list

iso_dic = {}
gene_list = []
outprots = open("%s/%s/%s_OGS_v2.1_prots.fasta" % (outdir, outprefix, outprefix), 'w')
reader = SeqIO.parse(inprots, format = 'fasta')
for rec in reader:
    if rec.id.split("-")[0] not in bac_genes:
        outprots.write(">%s\n%s\n" % (rec.id, str(rec.seq)))
        iso_dic[rec.id] = str(rec.seq)
        cur_gene = rec.id.split("-")[0]
        if cur_gene not in gene_list:
            gene_list.append(cur_gene)
outprots.close()

longest_list = longest_iso(gene_list, iso_dic)
outprots = open("%s/%s/%s_OGS_v2.1_prots_longest_isoform.fasta" % (outdir, outprefix, outprefix), 'w')
for longest in longest_list:
    outprots.write(">%s\n%s\n" % (longest, iso_dic[longest]))
outprots.close()

tran_dic = {}
outtrans = open("%s/%s/%s_OGS_v2.1_trans.fasta" % (outdir, outprefix, outprefix), 'w')
reader = SeqIO.parse(intranscripts, format = 'fasta')
for rec in reader:
    if rec.id.split("-")[0] not in bac_genes:
        outtrans.write(">%s\n%s\n" % (rec.id, str(rec.seq)))
        tran_dic[rec.id] = str(rec.seq)
outtrans.close()

outtrans = open("%s/%s/%s_OGS_v2.1_trans_longest_isoform.fasta" % (outdir, outprefix, outprefix), 'w')
for longest in longest_list:
    outtrans.write(">%s\n%s\n" % (longest, tran_dic[longest]))
outprots.close()
