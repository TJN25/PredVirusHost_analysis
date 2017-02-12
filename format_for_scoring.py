#!/usr/bin/python
'''This finds the columns containing the e values and then reassigns any that are less that 1e-9 to 
that threshold. The e values are then matched to the species/contigs in tab-delim format as arc	bac	euk	contig.
'''
input_file=open('output_with_species.txt','r') # the file that contains the species/contig, the protein and the e value
output_file=open('scores_calc.txt', 'a') 
output_file.write('Archaeal\tBacterial\tEukaryotic\tSpecies\n')

for line in input_file:
	words=line.rstrip()
	names=words.split('\t')
	arc=names[2][0:]
	bac=names[4][0:]
	euk=names[6][0:]
	try:
		arc=float(arc)
		bac=float(bac)
		euk=float(euk)
	except ValueError:
		print ('Error in line containing: %s %s %s' %(arc,bac,euk))
#	print arc
	if arc < 1e-9:
		arcx=1e-9
#		print arcx
	else:
		arcx=arc
	try:
		arcx=float(arcx)
	except ValueError:
		print ('Error in line containing: %s' % arcx)
		
#	print arcx
	if bac < 1e-9:
		bacx=1e-9
#		print bacx
	else:
		bacx=bac
	try:
		bacx=float(bacx)
	except:
		print ('Error in line containing: %s' % bacx)
#	print bacx
	if euk < 1e-9:
		eukx=1e-9
#		print eukx
	else:
		eukx=euk
	try:
		eukx=float(eukx)
	except ValueError:
		print ('Error in line containing: %s' % eukx)
#	print eukx
	if 'E-value' in words:
		pass
	else:
		try:
			species=names[7][0:]
			output_file.write('%s\t%s\t%s\t%s\n' % (arcx,bacx,eukx,species))	
		except IndexError:
			print words	
