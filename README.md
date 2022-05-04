# RNAG

RNAG - Block Gibbs Sampler for RNA Secondary Structure Prediction
version 1.2

Installation
Unpack the zip file. If you're reading this, you probably already have.
You will need some auxiliary programs in order to run RNAG:

probcons - available from http://probcons.stanford.edu/download.html<br />
Infernal package - available from http://eddylab.org/infernal/<br />
Vienna RNA package - available from http://www.tbi.univie.ac.at/RNA/<br />
Matlab or GNU Octave<br />
python 3.x<br />

Note: Matlab is a commercial product and requires a license.
Octave may be freely obtained from http://www.gnu.org/software/octave/
Octave is an free, open source alternative to Matlab. However,
it doesn't have all of the features of Matlab. In particular, Octave
seems lacking in the area of clustering. Currently, there are Octave
equivalents for pdist() and linkage() but not cluster() used in
hier_clus.m. A workaround for clustering is to use
mykmeans.m from http://www.christianherta.de/kmeans.html.

To use octave, edit rnag_init.py and change the line 'octave_prog = ...' to contain
the locatyion of teh octave executable. Change the line 'use_octave = 0' to 
'use_octave = 1'. 

Once the additional packages are installed, you will need to edit the file
rnag_init.py to point the variables to the correct locactions.

Running RNAG
To run, enter: 
python <path to rnag.py>
at the command prompt. This will display the command line options:

RNAG 1.2.0 2022-05-04
usage: python rnag.py fasta_file <iterations> <gamma_option> <Rfam_ID>
iterations - number of iterations, burn-in plus sample (optional, default = 1000)
gamma_option - 1 for range of gamma values (optional, default = 0)
Rfam_ID -  compare results to Rfam (optional, default = none)

Questions? Contact gibbs@brown.edu

Legalese:
Please acknowledge the program authors on any publication of scientific results
based in part on use of the program and cite the following article in which the
program was described.      

Wei, D., L. V. Alpert, et al. (2011). "RNAG: a new Gibbs sampler for predicting
RNA secondary structure for unaligned sequences." Bioinformatics 27(18):
2486-2493.

Copyright (C) 2011   
Brown University                                 
Providence, RI 02912
Email:  gibbs@brown.edu

RNAG is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.                                                                                 
RNAG is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
                                                                       
You should have received a copy of the GNU General Public License along with
 RNAG.  If not, see <http://www.gnu.org/licenses/>  or write to:
Free Software Foundation, Inc.
51 Franklin Street, Fifth Floor
Boston, MA  02110-1301, USA.                             
