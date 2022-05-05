#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
rnag.py - RNAG predict consensus secondary structures for unaligned sequences. 
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
import os.path
import re
import rnag_init
import fasta
import GS_seq

def print_options():
	print('usage: python rnag.py fasta_file <iterations> <gamma_option> <Rfam_ID>')
	print('iterations - number of iterations, burn-in plus sample (optional, default = 1000)')
	print('gamma_option - 1 for range of gamma values (optional, default = 0)')
	print('Rfam_ID -  compare results to Rfam (optional, default = none)')
	print()
	sys.exit(1)

def read_commandline():
	iter = 2000
	gamma_opt = 0
	rfam_id = None
	
	if len(sys.argv) == 1:
		print_options()
	elif len(sys.argv) == 2:
		fasta_file = sys.argv[1]
	elif len(sys.argv) == 3:
		fasta_file = sys.argv[1]
		iter = int(sys.argv[2])
	elif len(sys.argv) == 4:
		fasta_file = sys.argv[1]
		iter = int(sys.argv[2])
		gamma_opt = int(sys.argv[3])
	else: 
		fasta_file = sys.argv[1]
		iter = int(sys.argv[2])
		gamma_opt = int(sys.argv[3])
		rfam_id = sys.argv[4]

	if rfam_id is not None:
		p = re.compile('^RF\d{5}$')
		if not p.match(rfam_id):
			print('Error: ' + rfam_id + ' is not a valid RFam_Id')
			sys.exit(2)

	if iter % 2 == 1:
		iter += 1

	if iter < 6:
		print('Error: The number of iterations must be at least 6.')
		exit(3)
		
	return (fasta_file, iter, gamma_opt, rfam_id)


def xstr(s):
	if s is None:
		return 'none'
	return str(s)

def count_seqs(fasta_file):
	recs = fasta.read_fasta_file(fasta_file)
	return len(recs)


def initial_alignment(fasta_file):
	curr_dir = os.getcwd()
	if curr_dir[-1] != '/':
		curr_dir = curr_dir + '/'
	
	temp_file = curr_dir + 'aln.tmp'
	cmd = ' '.join([rnag_init.probcons_prog, fasta_file, '>', temp_file])
	os.system(cmd)

	recs = fasta.read_fasta_file(temp_file)
	
	align_file = curr_dir + '00aln'
	f = open(align_file, 'w')
	f.write('>alignment 0' + '\n')
	for rec in recs:
		f.write(rec.seq + '\n')
	f.close()

	os.remove(temp_file)
	if os.path.exists(curr_dir + '00str'):
		os.remove(curr_dir + '00str')


def read_rfam_file():
	cms = {}
	
	f = open(rnag_init.rfam_file)
	data = f.read()
	records = data.split('//\n')
	del records[-1]
	f.close()

	p = re.compile('ACCESSION (RF\d{5})')
	for rec in records:
		m = p.search(rec)
		acc = m.group(1)
		cms[acc] = rec

	return cms


def read_stockholm_file(file):
	struct = ''
	
	p = re.compile('^#=GC SS_cons')
	f = open(file)
	for line in f:
		if p.match(line):
			words = line.rstrip('\n').split()
			struct += words[2]

	p = re.compile('[<\[\{]')
	struct = p.sub('(', struct)
	p = re.compile('[>\]\}]')
	struct = p.sub(')', struct)
	p = re.compile('[^\(\)]')
	struct = p.sub('.', struct)

	return struct


def align_seqs(fasta_file, cm_file, curr_dir):
	structs = []
	seqs = fasta.read_fasta_file(fasta_file)
	for seq in seqs:
		seq_file = curr_dir + 'seq.tmp'
		f = open(seq_file, 'w')
		f.write('>' + seq.header + '\n')
		f.write(seq.seq + '\n')
		f.close()

		align_file = curr_dir + 'cm_align.tmp'
		temp_file = curr_dir + 'temp.tmp'

		cmd = ' '.join([rnag_init.cmalign_prog, '-o', align_file, cm_file, seq_file, '>', temp_file])
		os.system(cmd)
		struct = read_stockholm_file(align_file)
		structs.append(struct)
		
		os.remove(seq_file)
		os.remove(align_file)
		os.remove(temp_file)

	return structs

	
def setup_clustering(rfam_id, fasta_file):
	cms = read_rfam_file()
	if not rfam_id in cms:
		print('Error: ' + rfam_id + ' is not a valid RFam_Id')
		exit(3)

	curr_dir = os.getcwd()
	if curr_dir[-1] != '/':
		curr_dir = curr_dir + '/'
	cm_file = curr_dir + 'cm.tmp'
	f = open(cm_file, 'w')
	f.write(cms[rfam_id])
	f.write('//' + '\n')
	f.close()

	structs = align_seqs(fasta_file, cm_file, curr_dir)
	for i in range(len(structs)):
		struct_file = curr_dir + 'tru_str_' + str(i+1)
		f = open(struct_file, 'w')
		f.write(structs[i] + '\n')
		f.close()

	os.remove(cm_file)

##########################################################    

def main():
	VERSION = '1.2.0'
	DATE = '2022-05-04'
	print('RNAG ' + VERSION + ' ' + DATE)
	(fasta_file, iter, gamma_opt, rfam_id) = read_commandline()
	print('fasta file:   ' + fasta_file)
	print('iterations:   ' + str(iter))
	print('gamma option: ' + str(gamma_opt))
	print('Rfam ID:      ' + xstr(rfam_id))

	if rfam_id is not None:
		setup_clustering(rfam_id, fasta_file)
		
	initial_alignment(fasta_file)

	curr_dir = os.getcwd()
	if curr_dir[-1] != '/':
		curr_dir = curr_dir + '/'
	n_seqs = count_seqs(fasta_file)

	GS_seq.GS(curr_dir, fasta_file, n_seqs, range(1, n_seqs+1), iter)
	indicator_ref = 0
	if rfam_id is not None:
		indicator_ref = 1

	if rnag_init.rnag_path[-1] != '/':
		rnag_init.rnag_path = rnag_init.rnag_path + '/'

	if rnag_init.use_octave == 0:
		hier_clus_cmd = ' '.join(['hier_clus_web(',
								'\'' + curr_dir + '\'', ', ',
								str(n_seqs),  ', ',
								str(iter/2), ', ',
								str(iter/2), ', ',
								str(indicator_ref), ', ',
								str(gamma_opt),
								')'])
	else:
		hier_clus_cmd = ' '.join(['hier_clus_octave(',
								'\'' + curr_dir + '\'', ', ',
								str(n_seqs),  ', ',
								str(iter/2), ', ',
								str(iter/2), ', ',
								str(indicator_ref), ', ',
								str(gamma_opt),
								')'])   

	add_path_cmd = ''.join(["addpath('", 
							rnag_init.rnag_path,
							"');"])

	mat_file = curr_dir + 'matlab_cmd.m' 
	f = open(mat_file, 'w')
	f.write(add_path_cmd + '\n')
	f.write(hier_clus_cmd + '\n')
	f.close()

	if rnag_init.use_octave == 0:
		cmd = ' '.join([rnag_init.matlab_prog,  '-nosplash',
						'-nodisplay', '-nodesktop',
						'-r', "'matlab_cmd'"
						])
	else:
		cmd = ' '.join([rnag_init.octave_prog, 'matlab_cmd.m'])

	os.system(cmd)

if __name__ == "__main__":
	main()				
				
