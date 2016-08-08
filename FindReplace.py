#!/usr/bin/python
'''
Replaces the spaces in the input file with a single tab.
'''
def FindReplace(x,y)
input_file=open(x, 'r')
output_file=open(y, 'a')
for line in input_file:
	words=line.rstrip()
	words=words.replace(' - ',' ')
	names=words.replace(' ','\t')
	while '\t\t' in names:
		names=names.replace('\t\t','\t')
	output_file.write('%s\n'% names)
FindReplace('arVOG_file.txt','arVOG_file.tab')
FindReplace('baPOG_file.txt','baPOG_file.tab')
FindReplace('euVOG_file.txt','euVOG_file.tab')