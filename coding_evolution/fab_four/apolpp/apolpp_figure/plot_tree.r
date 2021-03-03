library(plotrix)
library(ape)

pdf("fig_apolpp_phylogeny.pdf", height = 5, width = 6.75)
par(mar = c(0,0,0,1), oma = c(0, 0, 0, 0), xpd = NA)

tree = read.tree("absrel_phylo.tree")
tree = drop.tip(tree, ">Frankliniella_occidentalis_XP_026284174.1") #this tip is misplaced easily -- out on a long branch
tree = ladderize(tree)
name_table = read.table("tip_names.txt", header = F, sep = "\t", stringsAsFactors=F, quote="")
tree$tip.label = name_table[[2]][match(tree$tip.label, name_table[[1]])]
tree$tip.label = sapply(tree$tip.label, function(x) parse(text=x))

soc_col = "#0099ff"

table = read.table("ape_edge_numbers.txt") #these are the branches that are detected by absrel -- translating between the numbers of nodes from absrel and the numbers of nodes from ape is a real pain in the ass
branch_cols = rep("black", length = 191)
for (edge in table$V1) {
    branch_cols[edge] = soc_col
}

plot(tree, edge.color = branch_cols, cex = 0.3, label.offset = 0.01)

text(2.55, 61.5, "Diptera", srt = 90)
segments(2.5, 56, 2.5, 67)

text(2.15, 74, "Lepidoptera", srt = 0)
segments(1.95, 69, 1.95, 79)

text(1.93, 84, "Coleoptera", srt = 0)
segments(1.74, 80, 1.74, 89)

text(1.78, 91.75, "Hemiptera", srt = 0)
segments(1.6, 90, 1.6, 94)

text(1.8, 10, "Halictidae", srt = 90)
segments(1.76, 1, 1.76, 19)

text(1.94, 15, "Bees", srt = 90)
segments(1.9, 1, 1.9, 29)

text(1.59, 38, "Ants", srt = 90)
segments(1.55, 30, 1.55, 47)

text(2.09, 30, "Hymenoptera", srt = 90)
segments(2.04, 1, 2.04, 55)

arrows(0.96, 16, 0.96, 19, length = 0.05, lwd = 2)

arrows(0.764, 32, 0.764, 35, length = 0.05, lwd = 2)

arrows(0.41, 53, 0.41, 50, length = 0.05, lwd = 2)

arrows(0.185, 64, 0.185, 67, length = 0.05, lwd = 2)

draw.circle(0.855, 22.4, 0.01)

dev.off()