## Jess: A 3D template searching program

Jess is a program that does 3D template searching on protein structures. The
idea is that the user provides one or more "templates", which consist of the 3D
coordinates of a group of atoms (e.g. a ligand binding motif), plus one or more
target structures in PDB format, and the program will try to find instances of each 
template on each target. The output is a PDB file constisting of the coordinates of
each hit, transformed so that they superpose onto the coordindates of the template,
and each hit is scored with a logE value, a probabilistic measure of goodness-of-fit.

For details on the algorithm, please read the [paper](10.1093/bioinformatics/btg226).  

Template nomenclature is based on the PDB format. For details on how to format templates,
please read the TESS (the predecesor of Jess) [paper](https://doi.org/10.1002%2Fpro.5560061104).

A detailed guide on template nomenclature will be released soon.

Last updated 16/4/2024

### Installation

Compile with something like:

`cd src`  
`gcc -c *.c`  
`gcc -o jess *.o -lm `  
`sudo mv jess /usr/local/bin`  

### Usage

`jess [template-list] [target-list] [rmsd] [distance] [max-dynamic-distance] [flags]`

* `template-list`: a list of filenames of TESS templates
* `target-list`: is a list of filenames of PDB files to search
* `rmsd`: the RMSD cutoff at which results are reported
* `distance`: the global distance cutoff used to guide the search
* `max-dynamic-distance`: maximum per-atom distance cutoff (details below). Set equal
                          to "distance" to override.

Important!
* As a rough estimate the distance cutoff should be 1.5-4 times
the RMSD cutoff. If you make it very large, the search becomes more
permissive, but execution time will suffer and the chance of finding spurious hits will increase.
* About the "max-dynamic-distance" argument: The user can define a dynamic matching distance cutoff
on a per-atom basis, that is added onto the global distance cutoff (if for example a single residue is flexible 
and is alowed to be matched with more relaxed spatial constraints). This dynamic distance of an atom can 
be optionally defined on the B-factor field of the ATOM record in the template. This argument is
the maximum allowed dynamic distance. To override dynamic distance completely, you can set this equal 
to the global distance argument.

[flags] : optional flags as a string with no spaces:  
* `f` : see PDB filenames in progress on stderr  
* `n` : do not transform coordinates of hit into	the template coordinate frame  
* `i` : include matches composed of residues belonging to
      multiple chains (if template is single-chain), or
	  matches with residues from a single chain
	  (if template has residues from multiple chains)  
* `q` : write filename of query instead of PDB ID from HEADER  
* `e` : parse atoms from all models separated by ENDMDL (use with
	  care). By default, Jess will only parse the first model

Example:

`cd examples`  
`../jess templates testfiles 2 3 3 > output`  

The output file is a flat file containing PDB fragments 
preceded by a single record of the form

`REMARK pdb-code rmsd template-file [debugging info]`

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

Each hit is followed by ENDMDL and a blank line.

### Filtering the output

Please note that in some cases, Jess performs multiple 
optimal aligments at a specific atom set in a given 
template-target pair. If you want to keep only the best hit,
pipe the output to the 'filter_jessout.py' script.

### Licence

The Jess code was written by Jonathan Barker ((c) 2002) 
and is now maintained and updated by Ioannis Riziotis [(e-mail)](mailto:ioannis.riziotis@crick.ac.uk).

This software was developed at EMBL-EBI in the [Thornton Group](https://www.ebi.ac.uk/research/thornton/)
and is provided for free under an [MIT License](https://choosealicense.com/licenses/mit/)
