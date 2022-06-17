•	README.TXT
WARNING: The custom scripts given in this folder are designed to work with the given LREC pipolins and thus might not work with other genbank pipolin files.
This folder contains the programs necessary for the generation of gggenomes input files, which can be easily obtained executing only the input_processing.py with Python 3 and the following requirements:
	1.	Software requirements:
		a.	Linux OS
		b.	Python3 with the following modules: Biopython, colormap
		c.	Minimap2
		d.	seq-gc from https://github.com/thackl/seq-scripts
	2.	File requirements
		a.	Genbank (.gbk) files in “original_gbk” folder
	b.	A file, order.txt, specifying the vertical order of the sequences when plotting them. Can be easily created with: ls original_gbk > order.txt (order by default).
The input_processing.py performs the following tasks:
	1.	Copies all gbk files from the original_gbk folder to the ordered_gbk folder including a prefix number that indicates the vertical order of the pipolin.
	2.	Modifies the gbk files in the ordered_gbk folder to make them valid for gggenomes
		a.	Changes the color of the features from RGB to HEX format
		b.	repeat_region features (atts) are given a unique identifier
		c.	Extract the plain sequence and create a new fasta file with it
	3.	Modified gbk and sequence files are stored in the processed_files folder.
	4.	Execute Synteny-GC_TIRs_script.bash, which performs the following tasks:
		a.	Calculate synteny with minimap2
		b.	Find TIRs (Terminal Inverted Repeats) with minimap2
		c.	Calculate CG content with seq-gc
After input_processing.py execution, the gggenomes input files generated are:
	•	renamed_files folder: contains the modified gbk and fasta files that will be used in gggenomes
	•	tirs.paf: paf files containing the terminal inverted repeats information and coordinates.
	•	gc.tsv: tsv file containing the gc content information.
	•	synteny.paf: paf file containing the synteny between sequences. 
Now, the gggenomes_pipolin_plotting.R script can be executed to obtain the figure. More information regarding the commands used are given inside the script. 
