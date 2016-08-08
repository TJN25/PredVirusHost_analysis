#!/usr/bin/python
'''
Formats the output of scores so that the 'rewrite_species_of_interest.py' can be used to make a 
file that can be graphed. The reason for that extra step is to prevent a graph for every sinlge 
contig being generated, when a lot of them will be blank. If you want all the contigs graphed then 
change line 10 from 'output_for_graph_metagenome.txt' to 'output_for_graph_species_of_interest.txt'.

'''
input_file=open('output_with_species.txt', 'r')
output_file=open('output_for_graph_metagenome.txt', 'a')
import math
for line in input_file:
	words=line.rstrip()
	names=words.split('\t')
	proteinname=names[0][0:]
	archaeal=names[2][0:]
	bacterial=names[4][0:]
	eukaryotic=names[6][0:]
	try:
		species=names[7][0:]
	except IndexError:
		pass
	archaeal=str(archaeal)
	archaeal=archaeal.replace('[','')
	archaeal=archaeal.replace(']','')
	archaeal=archaeal.replace("'","")
	archaeal=archaeal.replace(",","")
	archaeal=archaeal.replace(" ","\t")
	try:
		archaeal=float(archaeal)
		if archaeal==0:#		Needed because otherwise the log cannot be calculated.
			archaeal=350#		Assigned as this is slightly larger than the largest -logE that didnt = 0
		else:
			archaeal=math.log10(archaeal)
			archaeal=-archaeal
	except ValueError:
		pass
	bacterial=str(bacterial)
	bacterial=bacterial.replace('[','')
	bacterial=bacterial.replace(']','')
	bacterial=bacterial.replace("'","")
	bacterial=bacterial.replace(",","")
	bacterial=bacterial.replace(" ","\t")
	try:
		bacterial=float(bacterial)
		if bacterial==0:#		Needed because otherwise the log cannot be calculated.
			bacterial=350#		Assigned as this is slightly larger than the largest -logE that didnt = 0
		else:
			bacterial=math.log10(bacterial)
			bacterial=-bacterial		
	except ValueError:
		pass
	eukaryotic=str(eukaryotic)
	eukaryotic=eukaryotic.replace('[','')
	eukaryotic=eukaryotic.replace(']','')
	eukaryotic=eukaryotic.replace("'","")
	eukaryotic=eukaryotic.replace(",","")
	eukaryotic=eukaryotic.replace(" ","\t")
	try:
		eukaryotic=float(eukaryotic)
		if eukaryotic==0:#		Needed because otherwise the log cannot be calculated.
			eukaryotic=350#		Assigned as this is slightly larger than the largest -logE that didnt = 0
		else:
			eukaryotic=math.log10(eukaryotic)
			eukaryotic=-eukaryotic			
	except ValueError:
		pass
	if 'Protein Id' in proteinname:
		pass
	else:
		species=species.replace("'","")
		output_file.write('%s\t%s\tarchaeal\t%s\n%s\t%s\tbacterial\t%s\n%s\t%s\teukaryotic\t%s\n' % (proteinname,archaeal,species,proteinname,bacterial,species,proteinname,eukaryotic,species))
