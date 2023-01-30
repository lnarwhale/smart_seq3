#!/usr/bin/python
import string
import argparse

parser=argparse.ArgumentParser(description="cut of the smart-seq3-data")
parser.add_argument('-i','--input',help='the location of fastq before cut')
parser.add_argument('-o','--output',help='the location of output')
parser.add_argument('-n','--name',help='the name of fastq')
args=parser.parse_args()

input_file=args.input
output_file=args.output
fqname=args.name

gd1='CTAGTACGGGG'
gd2='TCGCCTTAGGG'
guo='CTGTCTCTTAT'
'''find umi and cut it'''
fil=open("{0}/deumi.fq".format(input_file),"r")
re1=open("{0}/sequ.fq".format(input_file),"w")
time=1
n1=0
n2=0
for pan in fil:
	loc1=pan.find(gd1)
	loc2=pan.find(gd2)
	if loc1 != -1:
		n1=n1+1
	if loc2 != -1:
		n2=n2+1
	time=time+1
	if time>20:
		break
if n1>n2:
	gd=gd1
else :
	gd=gd2

i=4
b=0
t=1
for li in fil:
	if i%4==0:
		name=li
	elif (i-1)%4==0:
		seq=li
		loc=seq.find(gd)
		if loc != -1:
			b=loc+len(gd)
			t=1
		else:
			b=0
			t=5
		if t <= 4:
			re1.write(name)
			re1.write(seq[b:])
	elif (i-2)%4==0:
		plu=li
		if t<= 4 :
			re1.write(plu)
	elif (i-3)%4==0:
		qua=li
		if t <= 4:
			re1.write(qua[b:])
	i=i+1
re1.close()
fil.close()

'''find CTGTCTCTTA and cut it'''
fil=open("{0}/sequ.fq".format(input_file),"r")
re1=open("{0}/re2.fq".format(input_file),"w")
i=4
b=0
for li in fil:
	if i%4==0:
		name=li
	elif (i-1)%4==0:
		seq=li
		loc=seq.find(guo)
		if loc != -1:
			b=loc-1
		else:
			b=len(seq)-1
		re1.write(name)
		re1.write(seq[:b]+'\n')
	elif (i-2)%4==0:
		plu=li
		re1.write(plu)
	elif (i-3)%4==0:
		qua=li
		re1.write(qua[:b]+'\n')
	i=i+1
re1.close()
fil.close()

'''to cut muti-G of 3'''
fil=open("{0}/re2.fq".format(input_file),"r")
res=open("{0}/{1}.fq".format(output_file,fqname),"w")
mark = 4
for line in fil:
	if mark%4==0 :
		name=line
		res.write(name)
	if (mark-1)%4==0 :
		seq=line
		n=2
		for i in seq:
			if n>len(seq):
				n=n-1
				break
			if seq[-n]=='G':
				n=n+1
			else:
				break
		res.write(seq[:-(n-1)]+'\n')
	if (mark-2)%4==0:
		plu=line
		res.write(plu)
	if (mark-3)%4==0:
		qua=line
		res.write(qua[:-(n-1)]+'\n')
	mark=mark+1
fil.close()
res.close()
