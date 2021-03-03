library(nlme)
library(geiger)
library(phytools)
library(ape)


behavs = read.table("behaviors.txt", row.names = 1, header = T)
table = na.omit(read.table("goterm_counts.txt", header = T, row.names = 1))
all_the_same = apply(table, 1, function(x) all(x == x[1]))
table = table[!all_the_same,]
table = table[rowSums(table) > 100,]
p_list = c()
goterm_list = c()
model_list = c()
for (i in 1:nrow(table)) {
    row = table[i,]
    goterm_list = c(goterm_list, row.names(row)[1])
    pgtab = merge(t(row), behavs, by = 0)
    colnames(pgtab) = c("taxon", "gocount", "behavior")
    rownames(pgtab) = pgtab$taxon
    tree = read.tree("RAxML_bestTree.halictid.tree")
    tree$edge.length = tree$edge.length * 100
    pglsModel = gls(gocount ~ behavior, correlation = corBrownian(phy = tree), data = pgtab, method = "ML")
    model_list = c(model_list, summary(pglsModel))
    cur_p = as.data.frame(anova(pglsModel))$p[2]
    p_list = c(p_list, cur_p)
}
corr_ps = p.adjust(p_list, method = "fdr")
outframe = cbind(goterm_list, signif(p_list,4), signif(corr_ps,4))
outframe = outframe[order(p_list),]
write.table(outframe, file = "alltaxa_pgls.txt", sep = "\t", quote = F, row.names = F, col.names = c("GOterm", "raw_p", "FDR_p"))



behavs = read.table("behaviors_simple.txt", row.names = 1, header = T)
table = na.omit(read.table("goterm_counts.txt", header = T, row.names = 1))
all_the_same = apply(table, 1, function(x) all(x == x[1]))
table = table[!all_the_same,]
table = table[rowSums(table) > 100,]
p_list = c()
goterm_list = c()
model_list = c()
for (i in 1:nrow(table)) {
  row = table[i,]
  goterm_list = c(goterm_list, row.names(row)[1])
  pgtab = merge(t(row), behavs, by = 0)
  colnames(pgtab) = c("taxon", "gocount", "behavior")
  rownames(pgtab) = pgtab$taxon
  tree = read.tree("RAxML_bestTree.halictid.tree")
  tree$edge.length = tree$edge.length * 100
  pglsModel = gls(gocount ~ behavior, correlation = corBrownian(phy = tree), data = pgtab, method = "ML")
  model_list = c(model_list, summary(pglsModel))
  cur_p = as.data.frame(anova(pglsModel))$p[2]
  p_list = c(p_list, cur_p)
}
corr_ps = p.adjust(p_list, method = "fdr")
outframe = cbind(goterm_list, signif(p_list,4), signif(corr_ps,4))
outframe = outframe[order(p_list),]
write.table(outframe, file = "alltaxa_binary_pgls.txt", sep = "\t", quote = F, row.names = F, col.names = c("GOterm", "raw_p", "FDR_p"))


behavs = read.table("behaviors_nopoly_noanc.txt", row.names = 1, header = T)
table = na.omit(read.table("goterm_counts.txt", header = T, row.names = 1))
table = subset(table, select=-c(HRUB,LALB,LCAL,MGEN,AVIR,NMEL,DNOV))
all_the_same = apply(table, 1, function(x) all(x == x[1]))
table = table[!all_the_same,]
table = table[rowSums(table) > 50,]
p_list = c()
goterm_list = c()
model_list = c()
for (i in 1:nrow(table)) {
  row = table[i,]
  goterm_list = c(goterm_list, row.names(row)[1])
  pgtab = merge(t(row), behavs, by = 0)
  colnames(pgtab) = c("taxon", "gocount", "behavior")
  rownames(pgtab) = pgtab$taxon
  tree = read.tree("RAxML_bestTree.halictid.tree")
  tree = drop.tip(tree, "MGEN")
  tree = drop.tip(tree, "DNOV")
  tree = drop.tip(tree, "HRUB")
  tree = drop.tip(tree, "LALB")
  tree = drop.tip(tree, "LCAL")
  tree = drop.tip(tree, "AVIR")
  tree = drop.tip(tree, "NMEL")
  tree$edge.length = tree$edge.length * 100
  pglsModel = gls(gocount ~ behavior, correlation = corBrownian(phy = tree), data = pgtab, method = "ML")
  model_list = c(model_list, summary(pglsModel))
  cur_p = as.data.frame(anova(pglsModel))$p[2]
  p_list = c(p_list, cur_p)
}
corr_ps = p.adjust(p_list, method = "fdr")
outframe = cbind(goterm_list, signif(p_list,4), signif(corr_ps,4))
outframe = outframe[order(p_list),]
write.table(outframe, file = "nopoly_noanc_pgls.txt", sep = "\t", quote = F, row.names = F, col.names = c("GOterm", "raw_p", "FDR_p"))




behavs = read.table("behaviors_nopoly_anc.txt", row.names = 1, header = T)
table = na.omit(read.table("goterm_counts.txt", header = T, row.names = 1))
table = subset(table, select=-c(HRUB,LALB,LCAL,MGEN))
all_the_same = apply(table, 1, function(x) all(x == x[1]))
table = table[!all_the_same,]
table = table[rowSums(table) > 50,]
p_list = c()
goterm_list = c()
model_list = c()
for (i in 1:nrow(table)) {
  row = table[i,]
  goterm_list = c(goterm_list, row.names(row)[1])
  pgtab = merge(t(row), behavs, by = 0)
  colnames(pgtab) = c("taxon", "gocount", "behavior")
  rownames(pgtab) = pgtab$taxon
  tree = read.tree("RAxML_bestTree.halictid.tree")
  tree = drop.tip(tree, "MGEN")
  tree = drop.tip(tree, "HRUB")
  tree = drop.tip(tree, "LALB")
  tree = drop.tip(tree, "LCAL")
  tree$edge.length = tree$edge.length * 100
  pglsModel = gls(gocount ~ behavior, correlation = corBrownian(phy = tree), data = pgtab, method = "ML")
  model_list = c(model_list, summary(pglsModel))
  cur_p = as.data.frame(anova(pglsModel))$p[2]
  p_list = c(p_list, cur_p)
}
corr_ps = p.adjust(p_list, method = "fdr")
outframe = cbind(goterm_list, signif(p_list,4), signif(corr_ps,4))
outframe = outframe[order(p_list),]
write.table(outframe, file = "nopoly_anc_pgls.txt", sep = "\t", quote = F, row.names = F, col.names = c("GOterm", "raw_p", "FDR_p"))

