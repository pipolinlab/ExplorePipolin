library(treeio)
library(ggtree)
library(dplyr)
library(ggplot2)

noShortCP017631 <- read.newick('~/Documents/tfm/the_whole_analysis/phylogeny/pi-polB/noShortCP017631/macse_noShortCP017631_NT_delexcl.fas.treefile')
ggtree(noShortCP017631) + geom_tiplab()

trimedAll <- read.newick('~/Documents/tfm/the_whole_analysis/phylogeny/pi-polB/all_trimmed/macse_trim_align_NT_delexcl.fas.treefile')
ggtree(trimedAll) + geom_tiplab()

p1 <- ggtree(noShortCP017631)
p2 <- ggtree(trimedAll)
d1 <- p1$data
d2 <- p2$data
d1$x <- d1$x
d2$x <- d2$x
d2$x <- d2$x + max(d1$x)
dd <- bind_rows(d1, d2) %>% filter(!is.na(label))
p1 + geom_tree(data = d2) + geom_tiplab(data = d2, size = 2) + 
  geom_line(aes(x, group = label, color = 'red'), 
            data = dd, alpha = 0.2, show.legend = F)
