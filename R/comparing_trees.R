# load libraries
library(phytools)
library(ggtree)
library(ggplot2)
#### read trees ####
pipolb_dir = '~/PycharmProjects/ExplorePipolin/data/new_analysis/phylogeny/piPolBs/'
core_dir = '~/PycharmProjects/ExplorePipolin/data/new_analysis/phylogeny/core/'
phy_pipolb_mafft <- read.newick(paste0(
  pipolb_dir, 'piPolB_concat.MAFFT-L-INS-i.renamed.txt.treefile'))
phy_pipolb_prank <- read.newick(paste0(
  pipolb_dir, 'piPolB_concat.PRANKcodon.renamed.txt.treefile'))
phy_core_mafft <- read.newick(paste0(
  core_dir, 'ROARYwMAFFT.core_gene_alignment.LREC_NCBI.aln.treefile'))

#### read and modify metadata ####
mdir <- '/home/liubov/PycharmProjects/ExplorePipolin/data/new_analysis/phylogeny/'
mdata <- read.csv(paste0(mdir, 'accessions_to_strainnames.tsv'),
                  sep = '\t', header = FALSE)
colnames(mdata) <- c('strain', 'accession', 'phylogroup', 'dataset')
mdata$color <- NA
mdata$color[mdata$phylogroup == 'A'] <- 'orange'
mdata$color[mdata$phylogroup == 'B1'] <- 'blue'
mdata$color[mdata$phylogroup == 'C'] <- 'green'
mdata$color[mdata$phylogroup == 'D'] <- 'sienna'


#### make a circular tree of core gene alignment ####
dataset_info = data.frame(strain = mdata$strain, dataset = mdata$dataset)

grp <- list(A = as.vector(subset(mdata$strain, mdata$phylogroup == 'A')),
            B1 = as.vector(subset(mdata$strain, mdata$phylogroup == 'B1')),
            C = as.vector(subset(mdata$strain, mdata$phylogroup == 'C')),
            D = as.vector(subset(mdata$strain, mdata$phylogroup == 'D')))

tree_core <- ggtree(phy_core_mafft, layout="circular", size = 0.5) + 
  geom_tiplab(size = 3.5, offset = 0.0005, show.legend = F) +
  geom_treescale(fontsize = 3, color = 'grey')

tree_core <- tree_core %<+% dataset_info + 
  geom_tippoint(aes(color = dataset), size = 1.5, show.legend = F)

groupOTU(tree_core, grp, 'Phylogroup') +
  aes(color = Phylogroup) +
  theme(legend.position = 'right') +
  scale_colour_manual(breaks = c('A','B1', 'C', 'D', 'NCBI', 'LREC'),
                      values=c('grey', 'orange', 'blue',
                               'green', 'sienna', 'red', 'black'))

#### phytools tanglegram ####
phy_tang <- cophylo(phy_core_mafft, phy_pipolb_mafft)

cols <- setNames(mdata$color, mdata$strain)

plot(phy_tang, fsize = 0.35, link.type='curved', link.lty = 'solid', link.lwd = 1,
     link.col = alpha(cols[phy_tang$assoc[,1]], 0.7), lwd = 1, pts = F, ftype = 'b')



#### dendextend tanglegram ####
library(DECIPHER)
library(dendextend)

dnd_pipolb_mafft <- ReadDendrogram(paste0(
  pipolb_dir, 'piPolB_concat.MAFFT-L-INS-i.renamed.txt.treefile'))
dnd_pipolb_prank <- ReadDendrogram(paste0(
  pipolb_dir, 'piPolB_concat.PRANKcodon.renamed.txt.treefile'))
dnd_core_mafft <- ReadDendrogram(paste0(
  core_dir, 'ROARYwMAFFT.core_gene_alignment.LREC_NCBI.aln.treefile'))

dndlist_pipolb_core <- dendlist(dnd_core_mafft, dnd_pipolb_mafft)
tanglegram(dndlist_pipolb_core, faster = FALSE, sort = TRUE, 
           margin_inner = 5, cex_main = 1.4,
           lab.cex = 0.8, lwd = 1,
           color_lines = cols[labels(dnd_core_mafft)],
           common_subtrees_color_branches = FALSE,
           highlight_branches_lwd = FALSE,
           main_left = 'Core genes', main_right = 'piPolBs',
           rank_branches = TRUE#, hang = TRUE
)

