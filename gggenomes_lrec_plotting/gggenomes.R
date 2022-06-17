#GGGenomes plotting script based on the instructions explained in https://thackl.github.io/gggenomes/index.html
#Make sure that the files (gbk, fa, paf, etc) are located in the correct directories
library(tidyverse)
library(gggenomes)
library(gtools)

### 1) Load GBK and FA files (features and sequences)
f <- read_feats(mixedsort(list.files("./processed_files/", "*.gbk$", full.names=TRUE)))
s <- read_seqs(mixedsort(list.files("./processed_files/", "*.fa$", full.names=TRUE)))

### 1.b) Simplest plot
p <- gggenomes(f, s) + geom_seq() + geom_gene()
p

### 2) Adding colors and more features
## Explanation: we can plot all the features together if we want to include the color code. 
## This is because the current (march 2022) version of gggenomes doesn't take pipolin specific features (repeat_regions, etc.)
## as genes, so the command "gggenomes(f, s) + geom_seq() + geom_gene(fill=f$hex_colour)" raises an error indicating
## that Aesthetics length (all features) is not the same as the data (only genes). 

##Therefore we subset the genes and later add the rest of features as additional tracks

### 2.A) Subset genes (CDS, rRNA, tRNA, tmRNA)
f_subset <- rbind(f[f$type=="CDS",], f[f$type=="rRNA",],
                 f[f$type=="tmRNA",], f[f$type=="tRNA",]) 
colour_vector <- paste("#", f_subset$colour, sep="") #Add "#" so the color code works properly

### 2.B) Add additional tracks (repeat_region, assembly_gap, etc.). 
#More custom tracks can easily be added copying a chunk of code and editing the feature type

#1st Subset features
#2nd Create a tibble with the seq_id, start and end fields. Color is not necessary as it is specified later
atts <- rbind(f[f$type=="repeat_region",])
atts_track <- tibble(seq_id=atts$seq_id,
                         start=atts$start, 
                         end=atts$end)

assembly_gap <- rbind(f[f$type=="assembly_gap",])
assembly_gap_track <- tibble(seq_id=assembly_gap$seq_id,
                             start=assembly_gap$start, 
                             end=assembly_gap$end)

misc_RNA <- rbind(f[f$type=="misc_RNA",])
misc_RNA_track <- tibble(seq_id=misc_RNA$seq_id,
                         start=misc_RNA$start, 
                         end=misc_RNA$end,
                         product=misc_RNA$product)

misc_feature <- rbind(f[f$type=="misc_feature",])
misc_feature_track <- tibble(seq_id=misc_feature$seq_id,
                         start=misc_feature$start, 
                         end=misc_feature$end)

partial_feature <- rbind(f[f$type=="partial-CDS",])
partial_feature_track <- tibble(seq_id=partial_feature$seq_id,
                             start=partial_feature$start, 
                             end=partial_feature$end)

## Note that some tracks may be empty because there is no such features among our pipolins. 
## We just don't include in the plot later. 

### 3) Include TIR, Synteny and GC content information

### 3.A) TIRs detection
tirs_paf <- read_paf("tirs.paf") %>% filter(seq_id == seq_id2 & start < start2 & map_length > 99 & de < 0.1)
tirs <- bind_rows(
  select(tirs_paf, seq_id=seq_id, start=start, end=end, de),
  select(tirs_paf, seq_id=seq_id2, start=start2, end=end2, de))

### 3.B) Synteny
firmicutes_links <- read_paf("synteny.paf")

### 3.C) CG content
firmicutes_gc <- thacklr::read_bed("gc.tsv")


### 4) Plot
## Plot has a ggplot-like syntax. We already saw the basic plot before. 
## Now we include the custom tracks, first inside feats=list(...), and later with +geom_feat(..)
## Empty tracks are not included (commented) to avoid errors

p <- gggenomes(f_subset, s, firmicutes_links, spacing = 0.5) + geom_seq() + geom_gene(fill=colour_vector)
p 

p <- (gggenomes(f_subset, s, feats=list(firmicutes_gc, assembly_gap_track, atts_track), #misc_RNA_track, tirs, partial_feature_track 
                firmicutes_links, spacing = 0.5)
      + geom_seq() + geom_gene(fill=colour_vector) 
      + geom_bin_label()
      + geom_feat(data=feats(atts_track), size=7, col="blue", alpha=0.7)
#      + geom_feat(data=feats(tirs), size=5, col="orange", alpha=0.5)
      + geom_feat(data=feats(assembly_gap_track), size=5, col="red", alpha=0.7)
#      + geom_feat(data=feats(misc_RNA_track), size=5, col="pink", alpha=1)
#      + geom_feat(data=feats(partial_feature_track), size=5, col="yellow", alpha=0.5)
      + geom_link(offset = c(0.3, 0.2))
#      + geom_gene_tag(aes(label=name), nudge_y=0.1, check_overlap = FALSE)
      )

p <- p + geom_ribbon(aes(x=(x+xend)/2, ymax=y+.24, ymin=y+.38-(.4*score),
                      group=seq_id, linetype="GC-content"), feats(firmicutes_gc),
                       fill="purple", position=position_nudge(y=-0.45))

p 

ggsave("LREC_EP_June_2022.svg", width=15, height=15, limitsize = FALSE)


