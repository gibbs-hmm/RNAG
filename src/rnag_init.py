#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
rnag_init.py - part of RNAG predict consensus secondary structures for unaligned sequences. 
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

# edit these to point to the proper paths

rna_fold_prog = '/mnt/d/BioComp/RNA/ViennaRNA_2_5/bin/RNAalifold'
cmbuild_prog = '/mnt/d/BioComp/infernal/infernal-1.0.2/bin/cmbuild'
cmalign_prog = '/mnt/d/BioComp/infernal/infernal-1.0.2/bin/cmalign'
probcons_prog = '/mnt/d/BioComp/RNA/RNAG/rnag_2022_05_04/probcons/probcons'
matlab_prog = '/home/bill/matlab_2022a/bin/matlab'
octave_prog = '/usr/bin/octave'  # octave location if octave is being used
rfam_file = '/mnt/d/BioComp/RNA/RNAG/rnag_2022_05_04/data/Rfam/Rfam.cm'
exit_cmd = '/mnt/d/BioComp/RNA/RNAG/rnag_2022_05_04/src/exit_cmd'  # for matlab
rnag_path = '/mnt/d/BioComp/RNA/RNAG/rnag_2022_05_04/src/'   # location where RNAG is installed

use_octave = 0   # change this to 1 to use octave instead of Matlab
