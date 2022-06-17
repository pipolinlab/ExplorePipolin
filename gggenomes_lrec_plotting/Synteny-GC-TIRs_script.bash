# 0) Create concatenated fasta
cat processed_files/*.fa > cat_seq.faa

# 1) Synteny
./minimap2/minimap2 -X -N 50 -p 0.1 -c cat_seq.faa cat_seq.faa > synteny.paf

# self-align opposite strands                        
for fna in `ls processed_files/*.fa`; do
  ./minimap2/minimap2 -c -B3 -O4 -E2 --rev-only $fna $fna > $fna.paf;
done;
cat processed_files/*paf > tirs.paf

#Empty ITRs output? no ITRs found

# 2) CG content
./seq-scripts-master/bin/seq-gc -Nbw 50 cat_seq.faa > gc.tsv
