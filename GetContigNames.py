#!/usr/bin/python
'''
Matches the protein and the contig name in a file that will be used to score each contig.
'''
import sys
input_file=open('fastafile.txt','r')
output_file=open('proteins_and_contigs_names.txt','a')
#delim=raw_input("Type the delimiter that is separating the contig name from the protein name (if unsure type 'help'): ")
delim=sys.argv[1]
if delim=='__':
	delim='|'
else:
	pass
#pos_start=raw_input("Type the number of delimiters that you must pass to be left with just the contig  (if unsure type 'help'): ")
pos_start=sys.argv[2]
#pos_end=raw_input("Type the number of delimiters from the end needed to be left with just the contig(if unsure type 'help'): ")
pos_end=sys.argv[3]
pos_start=int(pos_start)
pos_end=int(pos_end)
pos_end=-pos_end
for line in input_file:
	words=line.rstrip()
	if '>' in words:
		words=words.replace('>','')
		proteinname=words
		proteinname=proteinname.split(' ')
		proteinname=proteinname[0]
		proteinname=str(proteinname)
		proteinname=proteinname.replace("'",'')
		proteinname=proteinname.replace("[",'')
		proteinname=proteinname.replace(",",'')
		proteinname=proteinname.replace("]",'')
		proteinname=proteinname.replace(" ",'_')
		words=words.split(delim)
		if pos_end==0:
			contig=words[pos_start:][0:]
			contig=str(contig)
			contig=contig.replace("'",'')
			contig=contig.replace("[",'')
			contig=contig.replace(",",'')
			contig=contig.replace("]",'')		
			contig=contig.replace(" ",'_')		
			output_file.write('%s\t%s\n' %(proteinname,contig))
		else:
			contig=words[pos_start:pos_end][0:]
			contig=str(contig)
			contig=contig.replace("'",'')
			contig=contig.replace("[",'')
			contig=contig.replace(",",'')
			contig=contig.replace("]",'')		
			contig=contig.replace(" ",'_')		
			output_file.write('%s\t%s\n' %(proteinname,contig))			
	else:
		pass
