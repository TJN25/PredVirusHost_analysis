#!/usr/bin/python
input_file=open('fastafile.txt','r')
temp_file=open('temp_file.txt','a')

for line in input_file:
	line=line.rstrip()
	if '>' in line:
		words=line
		temp_file.write('\n%s ' % words)
	else:
		temp_file.write(line)
