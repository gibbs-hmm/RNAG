#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fasta.py - part of RNAG predict consensus secondary structures for unaligned sequences. 
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

class FastaRecord(object):
	def __init__(self, header, seq):
		self.header = header
		self.seq = seq


def read_fasta(infile):
	items = []
	count = 0;
	for line in infile:
		if line.startswith(">"):
			if count > 0:
				items.append(rec)
			header = line[1:].rstrip()
			seq = ""
			rec = FastaRecord(header,seq)
			count += 1
		else:
			seq += line.rstrip()
			rec = FastaRecord(header,seq)
	items.append(rec)
	
	return items


def read_fasta_file(filename):
	fd = open(filename)
	seqs = read_fasta(fd)
	fd.close()
	return seqs

