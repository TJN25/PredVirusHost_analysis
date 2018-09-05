#!/usr/bin/python
'''
Takes the output.txt file and removes any results that contain no significant hits.
'''
read_temp=open('output.txt','r')
output_file=open('output_sig_hits_only.txt','a')
for line in read_temp:
	words=line.rstrip()
	if 'No Significant Matches	1	No Significant Matches	1	No Significant Matches	1' in words:
		pass
	else:
		output_file.write('%s\n' % words)
#		print words		