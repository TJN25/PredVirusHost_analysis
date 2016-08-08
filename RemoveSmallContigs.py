#!/usr/bin/python
'''
remove_small_contigs.py looks at how many proteins are in each contig and removes any contigs from the fastafile if they
contain less than the number indicated by the user (default is 5 proteins or more).

'''
import sys
output_file=open('fastafile_subset.txt','a')
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
contig_cutoff=sys.argv[4]
contig_cutoff=int(contig_cutoff)
pos_start=int(pos_start)
pos_end=int(pos_end)
pos_end=-pos_end
dict_file=open('temp_file.txt','r')
x={}
for line in dict_file:
	line=line.rstrip()	
	protein=line.split()
	words=protein[0:1]
	words=str(words)
	words=words.split(delim)
	try:
		contig=words[pos_start:pos_end][0:]
		contig=str(contig)	
		contig=contig.replace("'",'')
		contig=contig.replace("[",'')
		contig=contig.replace(",",'')
		contig=contig.replace("]",'')		
		contig=contig.replace(" ",'_')	
		if contig in x:
			x[contig].append(protein)
		else:
			x[contig]=protein
	except IndexError:
		tmp_file=open('tmp1', 'a')
		words=str(words)
		tmp_file.write('%s\n' % words)
for key,item in x.items():
	key_length=len(item)
#	print key_length
	if key_length < contig_cutoff:
		del x[key]
	else:
		pass
for key,item in x.items():
	item=str(item)
	item=item.replace('[','')
	item=item.replace("'","")
	item=item.replace('"','')
	item=item.replace(',','')
	item=item.replace(']','')
	item=item.replace('>','\n>')
	item=item.replace('\n ','\n')
	item=item.replace(' ','\n')
	item=item.replace('\n\n','\n')
	output_file.write(item)
