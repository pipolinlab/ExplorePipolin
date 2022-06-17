
import os
from timeit import repeat
from colormap import rgb2hex
from Bio import SeqIO
import subprocess

if "processed_files" not in os.listdir():
	subprocess.run("mkdir processed_files", shell = True)
else:
	subprocess.run("rm -r processed_files", shell = True)
	subprocess.run("mkdir processed_files", shell = True)

if "ordered_gbk" not in os.listdir():
	subprocess.run("mkdir ordered_gbk", shell = True)
else:
	subprocess.run("rm -r ordered_gbk", shell = True)
	subprocess.run("mkdir ordered_gbk", shell = True)

### 0) Optional order
order_dict = {}
with open("order.txt", "r") as f:
    n_line = 0
    for line in f:
        n_line += 1
        order_dict[line.replace("\n","")] = str(n_line) 

for file in os.listdir("original_gbk"):
    if file in list(order_dict.keys()): order = order_dict[file]; subprocess.run("cp original_gbk/"+file+" ordered_gbk/"+order+"_"+file, shell=True)


### 1) Rename gbk and extract sequence
for file in os.listdir("ordered_gbk"): #change to input in final version
	num = file.split("_")[0]
	pipolin_id = file.split(".")[0].replace(num+"_","")
	repeat_n = 0
	with open("ordered_gbk/"+file, "r") as f:
		gbk_record = [record for record in SeqIO.parse(f, "genbank")][0]
		gbk_record.id = pipolin_id
		gbk_record.name = pipolin_id
		gbk_record.description = pipolin_id
		gbk_record.annotations['accessions'][0] = pipolin_id
		for feat in gbk_record.features:
			if 'colour' in list(feat.qualifiers.keys()):
				col_string = feat.qualifiers['colour'][0]
				col_list = col_string.split(" ")
				hex_col = rgb2hex(int(col_list[0]), int(col_list[1]), int(col_list[2]))
				feat.qualifiers['colour'] = [hex_col[1:]]
			if feat.type == 'repeat_region':
				repeat_n += 1
				feat.qualifiers['locus_tag'][0] = feat.qualifiers['locus_tag'][0]+"att"+str(repeat_n)
	with open("processed_files/processed_"+file, "w") as f:
		SeqIO.write(gbk_record, f, "genbank")
	
	### 2) Extract fasta from gbk
	fa_file_name = file.replace(".gbk", ".fa")
	SeqIO.convert("processed_files/processed_"+file, "genbank", "processed_files/processed_"+fa_file_name, "fasta")

### 3) Make concatenated fasta
subprocess.run("./Synteny-GC-TIRs_script.bash", shell=True)

