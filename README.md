Jess is Copyright (c) Jonathan Barker, 2002.
Now maintained and modified by Ioannis Riziotis (riziotis@ebi.ac.uk)

Last updated 14/12/2022

### INSTALLATION

Compile with something like:

`cd src`  
`gcc -c *.c`  
`gcc -o jess *.o -lm `  
`sudo mv jess /usr/local/bin`  

### OPERATING JESS

`jess [template-list] [target-list] [rmsd] [distance] [max-dynamic-distance] [flags]`

* `template-list`: a list of filenames of TESS templates
* `target-list`: is a list of filenames of PDB files to search
* `rmsd`: the RMSD cutoff at which results are reported
* `distance`: the global distance cutoff used to guide the search
* `max dynamic distance`: the maximum allowed template/query atom distance 
			  after adding the global distance cutoff and the 
			  individual atom distance cutoff defined in the
			  temperature field of the ATOM record in the template
			  file.

As a rough estimate the distance cutoff should be 1.5-4 times
the RMSD cutoff. But if you make it very large execution 
time will suffer. The smaller the better!

[flags]: optional flags as a string with no spaces:
	 f: see PDB filenames in progress on stderr
	 n: do not transform coordinates of hit into
	    the template coordinate frame
	 i: include matches composed of residues belonging to
	    multiple chains (if template is single-chain), or
	    matches with residues from a single chain
	    (if template has residues from multiple chains)
	 q: write filename of query instead of PDB ID from HEADER
	 e: parse atoms from all models separated by ENDMDL (use with
	     care). By default, Jess will only parse the first model

Example:

`cd examples`  
`../jess templates testfiles 2 3 3 > output`  

The output file is a flat file containing PDB fragments 
preceded by a single record of the form

`REMARK pdb-code rmsd template-file [some debugging info]`

where pdb-code is the PDB code of the hit file, rmsd is 
the rmsd after optimal superposition with the template and
template-file is the file containing the template which 
was hit.

The debugging info currently contains Det=number and
log(E)~number. If Det is not 1.0 then the superposition
is not valid (tell me about it please!). log(E) is a
preliminary statistical measure which should be used
with caution. Anything less than -4 should be a very
good hit and -3 is OK. (E is the expected number of hits
at random).

The PDB fragment which follows the remark is transformed 
into the coordinate frame of the template by optimal
superposition (if [n] flag is not there).

Each hit is followed by a blank line.

### FILTERING THE OUTPUT

Please note that in some cases, Jess performs multiple 
optimal aligments at a specific atom set in a given 
template-target pair. If you want to keep only the best hit,
pipe the output to the 'filter_jessout.py' scirpt.

### CAVEATS!

Jess and TESS rmsds are not the same!!!! This is because 
TESS aligns the site using a reference group, and Jess 
performs optimal superposition. In many cases this is 
desirable (see template_02 where there are many false
negatives in TESS, but few in Jess). The Jess rmsd will
*always* be <= the TESS rmsd. So cutoffs may well need
to be adjusted.

### BUGS!

Probably many. If you find one, tell me about it please!

