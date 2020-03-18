library(treeio)
library(ggtree)
library(dplyr)
library(ggplot2)
library(tidytree)
library(cowplot)
library(dendextend)
library(DECIPHER)

pipolb_dir = '~/Documents/tfm/the_whole_analysis/phylogeny/pi-polB/trying_tools/'
core_dir = '~/Documents/tfm/'
outgroups = c('18-1', 'MRY16-398')

#### roary_mafft vs roary_prank ####
core_mafft <- 
  read.newick(paste0(core_dir, 
                     'ROARYwMAFFT.core_gene_alignment.LREC_NCBI.aln.treefile'))
ggtree(core_mafft) + geom_tiplab(size = 1.5, show.legend = F)

core_prank <- 
  read.newick(paste0(core_dir,
                     'ROARYwPRANK.core_gene_alignment.LREC_NCBI.aln.treefile'))
ggtree(core_prank) + geom_tiplab(size = 1.5, show.legend = F)

# dendextend doesn't allow to create a tanglegram 
#  if the trees are not ultrametric
is.ultrametric(core_mafft)   # FALSE
is.ultrametric(core_prank)   # FALSE
# Not sure how this problem should be solved. I tried different ways 
#  to make the trees ultrametric. The following worked:
dnd_core_mafft <- 
  ReadDendrogram(paste0(core_dir,
                        'ROARYwMAFFT.core_gene_alignment.LREC_NCBI.aln.treefile'))
dnd_core_prank <- 
  ReadDendrogram(paste0(core_dir,
                        'ROARYwPRANK.core_gene_alignment.LREC_NCBI.aln.treefile'))
cor_cophenetic(dnd_core_mafft, dnd_core_prank)
## 0.9999
dndlist_cores <- dendlist(dnd_core_mafft, dnd_core_prank)
tanglegram(dndlist_cores, faster = FALSE, sort = TRUE, margin_inner = 5,
           lab.cex = 0.7, lwd = 1.5, #color_lines = c("blue"),
           common_subtrees_color_branches = TRUE,
           highlight_branches_lwd = FALSE
)

#### piPolB trees comparison ####
dnd_pipolb_mafft <- 
  ReadDendrogram(paste0(pipolb_dir,
                        'mafft/pi-polB_concat.MAFFT-L-INS-i.renamed.txt.treefile'))
dnd_pipolb_macse <- 
  ReadDendrogram(paste0(pipolb_dir,
                        'macse/macse_align_NT_delexcl.renamed.fas.treefile'))
dnd_pipolb_prank <- 
  ReadDendrogram(paste0(pipolb_dir,
                        'prank/pi-polB_concat.aln.best.renamed.fas.treefile'))

cor_cophenetic(dnd_pipolb_mafft, dnd_pipolb_prank, method_coef = "spearman")
## 0.849
cor_cophenetic(dnd_pipolb_macse, dnd_pipolb_mafft, method_coef = "spearman")
## 0.377
cor_cophenetic(dnd_pipolb_macse, dnd_pipolb_prank, method_coef = "spearman")
## 0.54

dndlist_mafft_prank <- dendlist(dnd_pipolb_mafft, dnd_pipolb_prank)
tanglegram(dndlist_mafft_prank, faster = FALSE, sort = TRUE, margin_inner = 5,
           lab.cex = 0.7, lwd = 1.5, #color_lines = c("blue"),
           common_subtrees_color_branches = TRUE,
           highlight_branches_lwd = FALSE
)
dndlist_macse_mafft <- dendlist(dnd_pipolb_macse, dnd_pipolb_mafft)
tanglegram(dndlist_macse_mafft, faster = FALSE, sort = TRUE, margin_inner = 5,
           lab.cex = 0.7, lwd = 1.5, #color_lines = c("blue"),
           common_subtrees_color_branches = TRUE,
           highlight_branches_lwd = FALSE
)
dndlist_macse_prank <- dendlist(dnd_pipolb_macse, dnd_pipolb_prank)
tanglegram(dndlist_macse_prank, faster = FALSE, sort = TRUE, margin_inner = 5,
           lab.cex = 0.7, lwd = 1.5, #color_lines = c("blue"),
           common_subtrees_color_branches = TRUE,
           highlight_branches_lwd = FALSE
)

#### core alignemnt vs piPolB aligments
cor_cophenetic(dnd_core_mafft, dnd_pipolb_mafft, method_coef = "spearman")
## 0.22
cor_cophenetic(dnd_core_mafft, dnd_pipolb_macse, method_coef = "spearman")
## 0.43
cor_cophenetic(dnd_core_mafft, dnd_pipolb_prank, method_coef = "spearman")
## 0.23

dndlist_core_pipolb_macse <- dendlist(dnd_core_mafft, dnd_pipolb_macse)
tanglegram(dndlist_core_pipolb_macse, faster = FALSE, sort = TRUE, margin_inner = 5,
           lab.cex = 0.7, lwd = 1.5, #color_lines = c("blue"),
           common_subtrees_color_branches = TRUE,
           highlight_branches_lwd = FALSE
)
dndlist_core_pipolb_mafft <- dendlist(dnd_core_mafft, dnd_pipolb_mafft)
tanglegram(dndlist_core_pipolb_mafft, faster = FALSE, sort = TRUE, margin_inner = 5,
           lab.cex = 0.7, lwd = 1.5, #color_lines = c("blue"),
           common_subtrees_color_branches = TRUE,
           highlight_branches_lwd = FALSE
)
dndlist_core_pipolb_prank <- dendlist(dnd_core_mafft, dnd_pipolb_prank)
tanglegram(dndlist_core_pipolb_prank, faster = FALSE, sort = TRUE, margin_inner = 5,
           lab.cex = 0.7, lwd = 1.5, #color_lines = c("blue"),
           common_subtrees_color_branches = TRUE,
           highlight_branches_lwd = FALSE
)

#### only reliable piPolB sequences ####
dnd_core_reliable <- 
  ReadDendrogram(paste0(core_dir,
                        'reliable.ROARYwMAFFT.core_gene_alignment.LREC_NCBI.aln.treefile'))
dnd_pipolb_reliable <- 
  ReadDendrogram(paste0(pipolb_dir,
                        'mafft/pi-polB_reliable.MAFFT-G-INS-i.renamed.txt.treefile'))
cor_cophenetic(dnd_core_reliable, dnd_pipolb_reliable, method_coef = "spearman")
## 0.208
dndlist_reliable_core_pipolb <- 
  dendlist(dnd_core_reliable, dnd_pipolb_reliable)
tanglegram(dndlist_reliable_core_pipolb, faster = FALSE, sort = TRUE, margin_inner = 5,
           lab.cex = 0.7, lwd = 1.5, #color_lines = c("blue"),
           common_subtrees_color_branches = TRUE,
           highlight_branches_lwd = FALSE
)

#### outgroups ####
core_OG_mafft <- 
  read.newick(paste0(core_dir,
                     'OGsp_ROARYwMAFFT.core_gene_alignment.LREC_NCBI_OG.aln.treefile'))
ggtree(core_OG_mafft) + geom_tiplab(size = 1.5,
                                    aes(color = label %in% outgroups),
                                    show.legend = F)
pipolb_OG_mafft <- 
  read.newick(paste0(pipolb_dir,
                     'mafft/OGsp_pi-polB_concat.MAFFT-L-INS-i.renamed.txt.treefile'))
ggtree(pipolb_OG_mafft) + geom_tiplab(size = 1.5,
                                      aes(color = label %in% outgroups),
                                      show.legend = F)

#### trimmed by macse ####
dnd_pipolb_macse_trimmed <- 
  ReadDendrogram(paste0(pipolb_dir,
                        'macse/trimmed/macse_trim_align_NT_delexcl.renamed.fas.treefile'))
dnd_pipolb_mafft_trimmed <- 
  ReadDendrogram(paste0(pipolb_dir,
                        'mafft/macse_trim_NT.MAFFT-L-INS-i.renamed.txt.treefile'))
dnd_pipolb_prank_trimmed <- 
  ReadDendrogram(paste0(pipolb_dir,
                        'prank/macse_trim_NT.aln.best.renamed.fas.treefile'))

cor_cophenetic(dnd_core_mafft, dnd_pipolb_macse_trimmed, method_coef = "spearman")
## 0.23
cor_cophenetic(dnd_core_mafft, dnd_pipolb_mafft_trimmed, method_coef = "spearman")
## 0.31
cor_cophenetic(dnd_core_mafft, dnd_pipolb_prank_trimmed, method_coef = "spearman")
## 0.22
# more or less the same results


# # ggtree
# p1 <- ggtree(core_mafft)
# p2 <- ggtree(core_prank)
# d1 <- p1$data 
# d2 <- p2$data
# d2$x <- d2$x + max(d1$x) + 0.0005
# dd <- bind_rows(d1, d2) %>% filter(!is.na(label))
# p1 + geom_tree(data = d2) + geom_tiplab(data = d2, size = 1.5) +
#   geom_line(aes(x, group = label, color = 'red'), 
#             data = dd, alpha = 0.2, show.legend = F)
# 
# p3 <- ggtree(pipolb_mafft)
# d3 <- p3$data
# d3$x <- d3$x + max(d1$x) + 0.0005
# ddd <- bind_rows(d1, d3) %>% filter(!is.na(label))
# p1 + geom_tree(data = d3) + geom_tiplab(data = d3, size = 1.5) +
#   geom_line(aes(x, group = label, color = 'red'), 
#             data = ddd, alpha = 0.2, show.legend = F)


# library(ape)
# comparePhylo(tree1, tree2, plot = TRUE, force.rooted = TRUE)
