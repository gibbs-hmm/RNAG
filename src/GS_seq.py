#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GS_seq.py - part of RNAG predict consensus secondary structures for unaligned sequences. 
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
import os
import parse_rna
import shutil
import rnag_init

def GS(Dir, fasta_file, num_seq, ind_seqs, iter, cons=0):    
	os.chdir(Dir)
	
	#A) Gibbs Sampling       
	#A.1) Initialization 1: from 00aln
	aa=open(Dir+'00aln','r')
	a=aa.readlines()
	aa.close()
	if os.path.exists('0.aln'):
		os.remove('0.aln')
		
	parse_rna.C_aln(Dir + str(0) + '.aln', fasta_file, a, num_seq)
	  
	for i in range(cons,iter):
		print(i)
		#1)RNAalifold:aln->str
		fold = os.popen(rnag_init.rna_fold_prog + ' -s 1 ' + Dir + str(i) + '.aln')
		#1.1)Record the consensus structure  - fixed problem with Vienna package 2.0
		result1=fold.readline()
		result2=fold.readline()
		result3=fold.readline()
		result = result3.split()
		if result3 == '':
			result = result2.split()
		result = result[0] + ' '
		#print result
		str_file=open(Dir+'00str','a')
		str_file.write('>structure '+str(i)+'\n')
#        str_file.write(result.strip()+'\n')
		str_file.write(result + '\n')
		str_file.close()
		#1.2)Prepare for Step 2)
		parse_rna.C_sto(Dir+str(i)+'.aln',Dir+str(i)+'.sto',result,num_seq,1)
		sto = os.popen(rnag_init.cmbuild_prog + ' -F '+ Dir + str(i) + '.cm '+ Dir + str(i) + '.sto')
		sys.stdout.flush()    
		for line in sto.readlines(): 
			print(line)
		#2)CM:str->aln
		cm = os.popen(rnag_init.cmalign_prog + ' --sample '+ Dir + str(i) + '.cm ' + fasta_file)
		ll=cm.readlines()
		#2.1)Record the consensus structure&Prepare for Step 1)
		alns=parse_rna.find_aln(Dir+str(i+1)+'.aln',ll,num_seq)
		aln_file=open(Dir+'00aln','a')
		aln_file.write('>>alignment  '+str(i+1)+'\n')
		for ii in range(0,num_seq):
			aln_file.write(alns[ii]+'\n')   
		aln_file.close()
		 #2.2)Clean up
		os.remove(Dir+str(i)+'.aln')        
		os.remove(Dir+str(i)+'.sto')
		if i!=iter-1:
			os.remove(Dir+str(i)+'.cm')        
		os.remove(Dir+'alirna.ps')
			
	#B) Matlab Analysis(which has nicer plots)
	#################################################################################
	for ind_seq in ind_seqs:
		parse_rna.project_strus(Dir+'project_'+str(ind_seq)+'.str',Dir+'00aln',Dir+'00str',ind_seq,num_seq) 

##iter=2000
### numbersequences = raw_input('Number of Sequences:')
##numbersequences = sys.argv[1]
##numbersequences = int(numbersequences)
##GS('./',numbersequences,range(1,numbersequences+1),iter)    
