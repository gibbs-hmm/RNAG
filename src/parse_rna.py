#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
parse_rna.py - part of RNAG predict consensus secondary structures for unaligned sequences. 
license: GPL 3

Copyright (C) 2022
Brown University                                 
Providence, RI 02912
Email:  gibbs@brown.edu

For details see:
Wei, D., L. V. Alpert, et al. (2011). "RNAG: a new Gibbs sampler for predicting
RNA secondary structure for unaligned sequences." Bioinformatics 27(18):
2486-2493.
"""

import sys

def chop_seq(seq,length=50):
	return [seq[i:i+length] for i in range(0, len(seq), length)]

def C_sto(read_aln,write_sto,stru,num_seq,spe):
	aln =open(read_aln,'r')
	out_pu =open(write_sto,'w')
	seqq=['']*num_seq
	out_pu.write("# STOCKHOLM 1.0\n\n");
	aln.readline()
	reader=aln.readline()
	while len(reader.strip())==0:
		reader=aln.readline()        
	if spe==0:
		strus=chop_seq(stru[0:stru.rfind(' ')],60)
	else:
		strus=chop_seq(stru[0:stru.rfind(' ')])
	chop_ind=0
	while reader!='' and reader[0]!=' ':
		for t in range(0,num_seq):
			sto_line=reader[0:reader.rfind(' ')]
			seqq[t]+=reader[30:reader.rfind(' ')]
			out_pu.write(sto_line+'\n')
			reader=aln.readline()     
		out_pu.write('#=GC SS_cons'+' '*18+strus[chop_ind]+'\n\n')
		chop_ind+=1
		reader=aln.readline()
		reader=aln.readline()                                  
	
	out_pu.write("//\n")
	out_pu.close()    
	aln.close()
	
def find_star(aln,num_seq):
	star=''
	for i in range(0,len(aln[0])):
		tmp=0
		for ii in range(num_seq-1):
			#print ii
			if aln[ii][i]!=aln[ii+1][i]:
				tmp+=1
		if tmp==0 and aln[0][i]!='-' and aln[0][i]!='.':
			star+='*'
		else:
			star+=' '
	return star
	
		   
	
def count(aln):
	tmp=aln.strip()
	cc=0
	for i in range(0,len(tmp)):
		if tmp[i]!='.' and tmp[i]!='-':
			cc+=1 
	return cc
	
def C_aln(output,in_seq,seqs,num_seq):
	seq_read=open(in_seq,'r')
	seq_names=seq_read.readlines()
	seq_read.close()
	aln_write=open(output,'w')
	aln_write.write('CLUSTAL 2.0.12 multiple sequence alignment\n\n')
	alns=[None]*num_seq
	sum=[0]*num_seq            
	for i in range(0,num_seq):
		alns[i]=chop_seq(seqs[i+1][0:-1])
	for i in range(0,len(alns[0])):    
		for ii in range(0,num_seq):
			name=seq_names[ii*2][1:-1]
			sum[ii]+=count(alns[ii][i])
#            print i, ii, name, '\n'
			aln_write.write(name+' '*(30-len(name))+alns[ii][i]+' '+str(sum[ii])+'\n')
		aln_write.write(' '*30+find_star([row[i] for row in alns],num_seq)+'\n\n')               
							
	aln_write.close()
	
	
def find_aln(to,stru,num_seq):
	i=0
	alns=num_seq*['']
	aln_write=open(to,'w')
	aln_write.write('CLUSTAL 2.0.12 multiple sequence alignment\n\n')        
	seqs=['']*num_seq
	len_seq=[0]*num_seq
	while '# STOCKHOLM' not in stru[i]:
		i+=1
		if i >= len(stru):
			f = open('error.txt', 'w')
			f.write('Error: cmalign count not build alignment\n')
			f.close()
			sys.exit('Error: cmalign count not build alignment')
	else:
		i+=3
		while '#\n' not in stru[i]:
			for ii in range(0,num_seq):
				tmp=stru[i]
				seqs[ii]=tmp[tmp.rfind(' ')+1:-1]
				len_seq[ii]+=count(seqs[ii])
				aln_write.write(tmp[0:tmp.rfind(' ')]+' '*(30-tmp.rfind(' '))+tmp[tmp.rfind(' ')+1:-1]+' '+str(len_seq[ii])+'\n')                                       
				alns[ii]+=seqs[ii]
				i+=1
			i+=3
			aln_write.write(' '*30+find_star(seqs,num_seq)+'\n\n')
	aln_write.close() 
	return alns                

def project(in_aln,in_str):
	tt=[1]*len(in_aln)
	output=in_str
	indent=[]
	for ii in range(0,len(in_aln)):
		if in_aln[ii]=='-' or in_aln[ii]=='.':
			indent.append(ii)
			tt[ii]='0'
	while indent:
		#print indent[0]
		if output[indent[0]]!='.':
			#print indent[0]
			if in_str[indent[0]]=='(' or in_str[indent[0]]=='<':
				left=1
				right=0
				step=1
			else:
				left=0
				right=1   
				step=-1
			count=indent[0]
			output=output[:count] +'.'+ output[count+1:]
			#print output
			while left!=right:
				count=count+step
				if in_str[count]=='(' or in_str[count]=='<':
					left=left+1
				elif in_str[count]==')' or in_str[count]=='>':
					right=right+1;
			output=output[:count] +'.'+ output[count+1:]
			#print output
		del indent[0]
	out2=''
	for ww in range(len(output)):
		if tt[ww]==1:
			out2+=output[ww]
	return out2    

def project_strus(to,input_aln,input_str,ind_seq,num_seq):
	stru_write=open(to,'w')
	stru_read=open(input_str)
	stru=stru_read.readlines()
	stru_read.close()
	aln_read=open(input_aln)
	aln=aln_read.readlines()
	aln_read.close()
	N=int(len(stru)/2)
	list=[[]]*N
	#print N
	for i in range(N):
		print(i)
		out2=project(aln[(num_seq+1)*i+ind_seq][0:-1],stru[1+2*i].strip())        
		list[i]=bra2list(out2)
		stru_write.write('>Structrue '+str(i)+'\n')
		stru_write.write(out2+'\n')            
	stru_write.close()
	return list

def bra2list(bracket):
	N=len(bracket)
	list=[]
	topstack = -1
	openstack = []    
	for i in range(N):
		if bracket[i]=='(':
			topstack+=1
			openstack.append(i) 
		elif bracket[i]==')':
			list.append(openstack[topstack]*N+i) 
			del openstack[topstack]
			topstack = topstack-1
	return list
					
					
def bratrans(tt):
	truth=''
	if(tt[len(tt)-1]=='\n'):
		tt=tt[:-1]
	for i in range(len(tt)):
		if(tt[i]=='<'):
			truth+='('
		elif(tt[i]=='>'):
			truth+=')'
		else:
			truth+='.'
	return truth

def isthere(listt,num):
	yes=0
	for i in listt:
		if(i==num):
			yes=1
			break
	return yes

def sta(stru,true,N=1):
	TP=0
	for i in range(len(stru)):
		TP+=isthere(true,stru[i])
	print(TP,len(stru),len(true))
	TPR=float(TP)/(float(len(true)))
	if len(stru)==0:
		PPV=0
	else:
		PPV=float(TP)/(float(len(stru)))    
	#FP=len(stru)-TP
	#FN=len(true)-TP;
	#TN=N*(N-1)/2-len(true)-FP
	#TPR=float(TP)/(float(TP+FN))
	#PPV=float(TP)/(float(TP+FP))
	#FPR=FP/(FP+TN)
	sta=[PPV,TPR]
	return sta
	

