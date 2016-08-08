#!/usr/bin/python
'''
This allows the results from contigs that you want, to be graphed in R. There is a sample output avaialable.

You will need to make a list of the species you want to graph, including the extra '[' and other
things found in the output_with_species.txt file and save into the same directory and name it 
'species_to_graph.txt'
'''
input_file=open('species_to_graph.txt', 'r')
input_file_2=open('output_for_graph_metagenome.txt', 'r')
output_file=open('output_for_graph_species_of_interest.txt', 'a')
a={}
for line in input_file:
	line=line.rstrip()
	value=1
	a[line]=value
for line in input_file_2:
	line=line.rstrip()
	words=line.split('\t')
	names=words[3]
	if names in a:
		output_file.write('%s\n' % line)
	else:
		pass
