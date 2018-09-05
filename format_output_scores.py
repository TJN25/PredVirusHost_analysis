#!/usr/bin/python
'''
Reformats the output from the r script into a tab delimited file that is easily searchable.
'''
input_file=open('output_scores.txt','r')
output_file=open('scores.txt','a')
output_file.write('Contig\tArchaeal_P_Value\tBacterial_P_Value\tEukaryotic_P_Value\tArchaeal_Mean\tBacterial_Mean\tEukaryotic_Mean\tArchaeal_Median\tBacterial_Median\tEukaryotic_Median\tNumber_of_Proteins')
for line in input_file:
	line=line.rstrip()
	names=line.replace(' ','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	names=names.replace('\t\t','\t')
	words=names	
	words=words.replace('[','\n[')
	output_file.write(words)
	
	
