#!/usr/bin/python
'''
Adds the name of the contig that each protein belongs to into a new column in the output.txt file.
'''
dict_file=open('proteins_and_contigs_names.txt','r')
input_file=open('output.txt', 'r')
output_file=open('output_with_species.txt','a')
output_file.write("'Protein Id'\t 'Archaeal Virus Model'\tE-value\t'Bacterial Phage Model'\tE-value\t'Eukaryotic Virus Model'\tE-value\tContig\n")
x={}
for line in dict_file:
	words=line.rstrip()
	names=words.split('\t')
	key=names[0][0:]
	value=names[1:][0:]
	key=str(key)
	item=key.replace('[','')
	item=item.replace("'","")
	item=item.replace(',','')
	key=item.replace(']','')	
	x[key]=value
for line in input_file:
	words=line.rstrip()
	names=words.split('\t')
	key=names[0][0:]
	key=str(key)
	item=key.replace('[','')
	item=item.replace("'","")
	item=item.replace(',','')
	key=item.replace(']','')		
	if key in x:
		species=x[key]	
		output_file.write('%s\t%s\n' % (words,species))
	else:
		print('%s output_with_contigs.py' % key)