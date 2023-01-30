#!/usr/bin/python
import string
import argparse

parser=argparse.ArgumentParser(description="precut of the smart-seq3-data")
parser.add_argument('-i','--input',help='the location of fastq before cut')
parser.add_argument('-o','--output',help='the location of output')
parser.add_argument('-n','--name',help='the name of fastq')
args=parser.parse_args()

input_file=args.input
output_file=args.output
name=args.name

um1='CTAGTACGGGG'
um2='TCGCCTTAGGG'
guo='CTGTCTCTTAT'
'''find umi and cut it'''
fil=open("{0}/{1}.fq".format(input_file,name),"r")
re1=open("{0}/umi.fq".format(output_file),"w")
time=1
n1=0
n2=0
for pan in fil:
    loc1=pan.find(um1)
    loc2=pan.find(um2)
    if loc1 != -1:
        n1=n1+1
    if loc2 != -1:
        n2=n2+1
    time=time+1
    if time>20:
        break
if n1>n2:
    umi=um1
else :
    umi=um2

i=4
b=0
t=1
for li in fil:
    if i%4==0:
        name=li
    elif (i-1)%4==0:
        seq=li
        loc=seq.find(umi)
        if loc != -1:
            b=loc-8
            if b<0:
                t=5
            else:	 
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


