#!/usr/bin/python
'''
Takes the output of the hmm search and changes it into a format that compares the e-values of each protein match to each superkingdom. 
A full output file containing all proteins.
'''


import math
ar_input=open('arVOG_list.txt','r')
ba_input=open('baPOG_list.txt','r')
eu_input=open('euVOG_list.txt','r')
protein_list=open('proteinnames_list.txt','r')
temp_file=open('output_tmp.txt', 'a')
ar={}
for line in ar_input:
	words=line.rstrip()
	names=words.split()
	key=names[0][0:]
#	print key
	model=names[1][0:]
	arvalue=names[2][0:]
	arvalue=arvalue.split()
	y=arvalue
#	print model
#	print value
	
	if key in ar:
		x=ar[key]
		try:
			x=x[1][0:]
			x=str(x)
			x=x.replace('[','')
			x=x.replace(']','')
			x=x.replace("'","")
			x=x.replace(",","")
			x=x.replace(" ","\t")
			try:
				x=float(x)
				if x==0:#		Needed because otherwise the log cannot be calculated.
					xvalue=350#		Assigned as this is slightly larger than the largest -logE that didnt = 0
				else:
					x=math.log10(x)
					xvalue=-x
			except ValueError:
				print ('Ignoring line containing: %s' % x)
			y=str(y)
			y=y.replace('[','')
			y=y.replace(']','')
			y=y.replace("'","")
			y=y.replace(",","")
			y=y.replace(" ","\t")
			try:
				y=float(y)
				if y==0:#		Needed because otherwise the log cannot be calculated.
					yvalue=350#		Assigned as this is slightly larger than the largest -logE that didnt = 0
				else:
					y=math.log10(y)
					yvalue=-y
			except ValueError:
				print ('Ignoring line containing: %s' % y)
			try:
				if xvalue > yvalue:
					pass
				else:
					del ar[key]
					model=model.split()
					ar[key]=model
#					value=value.split()
					ar[key].append(arvalue)
			except NameError:
				print 'No matches'
		except IndexError:
			pass
	else:
		model=model.split()
		ar[key]=model
		ar[key].append(arvalue)
ba={}
for line in ba_input:
	words=line.rstrip()
	names=words.split()
	key=names[0][0:]
#	print key
	model=names[1][0:]
	model=model.replace('bact','ba')
	bavalue=names[2][0:]
	bavalue=bavalue.split()
	y=bavalue
#	print model
#	print value
	
	if key in ba:
		x=ba[key]
		try:
			x=x[1][0:]
			x=str(x)
			x=x.replace('[','')
			x=x.replace(']','')
			x=x.replace("'","")
			x=x.replace(",","")
			x=x.replace(" ","\t")
			try:
				x=float(x)
				if x==0:#		Needed because otherwise the log cannot be calculated.
					xvalue=350#		Assigned as this is slightly larger than the largest -logE that didnt = 0
				else:
					x=math.log10(x)
					xvalue=-x
			except ValueError:
				print ('Ignoring line containing: %s' % x)
			y=str(y)
			y=y.replace('[','')
			y=y.replace(']','')
			y=y.replace("'","")
			y=y.replace(",","")
			y=y.replace(" ","\t")
			try:	
				y=float(y)
				if y==0:#		Needed because otherwise the log cannot be calculated.
					yvalue=350#		Assigned as this is slightly larger than the largest -logE that didnt = 0
				else:
					y=math.log10(y)
					yvalue=-y
			except ValueError:
				print ('Ignoring line containing: %s' % y)
			try:
				if xvalue > yvalue:
					pass
				else:
					del ba[key]
					model=model.split()
					ba[key]=model
#					value=value.split()
					ba[key].append(bavalue)
			except NameError:
				print 'No matches'
		except IndexError:
			pass
	else:
		model=model.split()
		ba[key]=model	
		ba[key].append(bavalue)
eu={}
for line in eu_input:
	words=line.rstrip()
	names=words.split()
	key=names[0][0:]
#	print key
	model=names[1][0:]
	euvalue=names[2][0:]
	euvalue=euvalue.split()
	y=euvalue
#	print model
#	print value
	
	if key in eu:
		x=eu[key]
		try:
			x=x[1][0:]
			x=str(x)
			x=x.replace('[','')
			x=x.replace(']','')
			x=x.replace("'","")
			x=x.replace(",","")
			x=x.replace(" ","\t")
			try:
				x=float(x)
				if x==0:#		Needed because otherwise the log cannot be calculated.
					xvalue=350#		Assigned as this is slightly larger than the largest -logE that didnt = 0
				else:
					x=math.log10(x)
					xvalue=-x
			except ValueError:
				print ('Ignoring line containing: %s' % x)
			y=str(y)
			y=y.replace('[','')
			y=y.replace(']','')
			y=y.replace("'","")
			y=y.replace(",","")
			y=y.replace(" ","\t")
			try:
				y=float(y)
				if y==0:#		Needed because otherwise the log cannot be calculated.
					yvalue=350#		Assigned as this is slightly larger than the largest -logE that didnt = 0
				else:
					y=math.log10(y)
					yvalue=-y
			except ValueError:
				print ('Ignoring line containing: %s' % y)
			try:
				if xvalue > yvalue:
					pass
				else:
					del eu[key]
					model=model.split()
					eu[key]=model
#					value=value.split()
					eu[key].append(euvalue)
			except NameError:
				print 'No matches'
		except IndexError:
			pass
	else:
		model=model.split()
		eu[key]=model
		eu[key].append(euvalue)
for line in protein_list:
	words=line.rstrip()
	if words in ar:
		arvalue=ar[words]
		arvalue=str(arvalue)
		arvalue=arvalue.replace('[','')
		arvalue=arvalue.replace(']','')
		arvalue=arvalue.replace("'","")
		arvalue=arvalue.replace(",","")
		arvalue=arvalue.replace(" ","\t")	
	else:
		arvalue="'No Significant Matches'\t1"
	if words in eu:
		euvalue=eu[words]
		euvalue=str(euvalue)
		euvalue=euvalue.replace('[','')
		euvalue=euvalue.replace(']','')
		euvalue=euvalue.replace("'","")
		euvalue=euvalue.replace(",","")
		euvalue=euvalue.replace(" ","\t")
	else:
		euvalue="'No Significant Matches'\t1"
	if words in ba:
		bavalue=ba[words]
		bavalue=str(bavalue)
		bavalue=bavalue.replace('[','')
		bavalue=bavalue.replace(']','')
		bavalue=bavalue.replace("'","")
		bavalue=bavalue.replace(",","")
		bavalue=bavalue.replace(" ","\t")
	else:
		bavalue="'No Significant Matches'\t1"
	temp_file.write('%s\t%s\t%s\t%s\n' % (words,arvalue,bavalue,euvalue))

