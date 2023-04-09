import utils
from Gene import Gene
import cPickle as pickle

pickle_file = open("/scratch/tmp/berubin/resequencing/hic/LOEN/genotyping/LOEN_gene_vcf_dic.pickle", 'rb')
gene_objects = pickle.load(pickle_file)
pickle_file.close()

loen_genes = ["02440", "07629", "10519", "03154", "05010", "05308", "09405", "10614"]
lpau_genes = ["05479", "08784", "10601", "06573", "11655", "16761", "03601", "10928"]
new_dic = {}
for gene in loen_genes:
    new_dic["LOEN_%s" % gene] = gene_objects["LOEN_%s" % gene]
new_dic["LOEN_10709"] = gene_objects["LOEN_10709"]
new_dic["LOEN_07008"] = gene_objects["LOEN_07008"]
pickle_file = open("/scratch/tmp/berubin/resequencing/hic/LOEN/genotyping/testset_LOEN_gene_vcf_dic.pickle", 'wb')
pickle.dump(new_dic, pickle_file)
pickle_file.close()

pickle_file = open("/scratch/tmp/berubin/resequencing/hic/LOEN/genotyping/testset_LOEN_gene_vcf_dic.pickle", 'rb')
gene_objects = pickle.load(pickle_file)
print gene_objects
for k, g in gene_objects.items():
    print g.name
    print g.scaf
    print g.fourfold


pickle_file = open("/scratch/tmp/berubin/resequencing/hic/LPAU/genotyping/LPAU_gene_vcf_dic.pickle", 'rb')
gene_objects = pickle.load(pickle_file)
pickle_file.close()

new_dic = {}
for gene in lpau_genes:
    new_dic["LPAU_%s" % gene] = gene_objects["LPAU_%s" % gene]

new_dic["LPAU_05002"] = gene_objects["LPAU_05002"]
new_dic["LPAU_01477"] = gene_objects["LPAU_01477"]
pickle_file = open("/scratch/tmp/berubin/resequencing/hic/LPAU/genotyping/testset_LPAU_gene_vcf_dic.pickle", 'wb')
pickle.dump(new_dic, pickle_file)
pickle_file.close()

pickle_file = open("/scratch/tmp/berubin/resequencing/hic/LPAU/genotyping/testset_LPAU_gene_vcf_dic.pickle", 'rb')
gene_objects = pickle.load(pickle_file)
print gene_objects
for k, g in gene_objects.items():
    print g.name
    print g.scaf
    print g.fourfold
