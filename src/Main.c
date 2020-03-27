// ==================================================================
// Main.c
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// A simple controlling shell for Jess
// ==================================================================

#include "Jess.h"
#include "TessTemplate.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>

// ==================================================================
// Global constants
// ==================================================================
// atomFormat			The format of a PDB ATOM record (for printf)
// ==================================================================

static const char *atomFormat =
	"ATOM  %5i%5s%c%-3s%c%c%4i%-4c%8.3f%8.3f%8.3f\n"; //Riziotis edit
	//"ATOM  %5i%5s%c%-4s%c%4i%-4c%8.3f%8.3f%8.3f\n";

// ==================================================================
// Global flags
// ==================================================================
// feedbackQ			Give feedback while processing
// ==================================================================

static int feedbackQ=0;

// ==================================================================
// Local functions
// ==================================================================

static void output(
	const Atom *A,
	const double *M,
	const double *c,
	const double *v,
	int no_transform
	)
{
	double x[3];
	int i,j,k;
	char name[5];
	char resName[4];

	// First transform the coordinates according to the
	// transform and centroids given as arguments. See
	// also Super.h/Super.c (type Superposition).

	//Riziotis edit
	if (no_transform == 1)
	{
		for(i=0; i<3; i++)
		{
			x[i]=A->x[i];
		}
	}
	
	else
	{
		for(i=0; i<3; i++)
		{
			x[i]=0.0;
			for(j=0; j<3; j++)
			{
				x[i] += M[3*i+j]*(A->x[j]-c[j]);
			}

			x[i] += v[i];
		}
	}

	// Remove underscores from atom name and residue name

	strncpy(name,A->name,4);
	strncpy(resName,A->resName,3);
	name[4]=0;
	resName[3]=0;
	for(i=0; i<3; i++)
	{
		if(resName[i]=='_') resName[i]=' ';
	}
	for(i=0; i<4; i++)
	{
		if(name[i]=='_') name[i]=' ';
	}

	// Output the ATOM record with the coordinates
	// and names suitably transformed.

	printf(
		atomFormat,
		A->serial,
		name,
		A->altLoc,
		resName,
		//Riziotis edit
		A->chainID1,
		A->chainID2,
		//A->chainID,
		A->resSeq,
		A->iCode,
		x[0],x[1],x[2],
		A->occupancy,
		A->tempFactor,
		A->segID,
		A->element,
		A->charge
		);
}
static void search(const char *filename,Jess *J,double tRmsd,double tDistance, int no_transform)
{
	Molecule *M;
	Superposition *sup;
	Template *T;
	Atom **A;
	FILE *file;
	JessQuery *Q;
	int i,j,k,count;
	const double *P,*c[2];
	double det;
	double logE;

	if(!(file=fopen(filename,"r")))
	{
		perror(filename);
		return;
	}

	M = Molecule_create(file);

	fclose(file);
	if(!M)
	{
		fprintf(stderr,"%s: bad PDB file\n",filename);
		return;
	}

	Q=Jess_query(J,M,tDistance);

	while(JessQuery_next(Q))
	{
		T=JessQuery_template(Q);

		count=T->count(T);
		sup = JessQuery_superposition(Q);
		A = JessQuery_atoms(Q);

		if(Superposition_rmsd(sup)<=tRmsd)
		{
			P=Superposition_rotation(sup);

			// This is to check for the unusual case where
			// P is computed as a matrix with determinant -1.
			// THIS SHOULD NOT HAPPEN!

			det = 0.0;
			det += P[0]*(P[4]*P[8]-P[5]*P[7]);
			det -= P[1]*(P[3]*P[8]-P[5]*P[6]);
			det += P[2]*(P[3]*P[7]-P[4]*P[6]);

			c[0]=Superposition_centroid(sup,0);
			c[1]=Superposition_centroid(sup,1);


			logE=T->logE(T,Superposition_rmsd(sup),Molecule_count(M));

			printf("REMARK %s ",Molecule_id(M) ? Molecule_id(M):filename);
			printf("%.3f ",Superposition_rmsd(sup));
			printf("%s Det= %.1f log(E)~ %.2f\n",T->name(T),det,logE);

			// Output the transformed target atoms if reverseQ is
			// not specified.

			for(i=0; i<count; i++)
			{
				output(A[i],P,c[0],c[1],no_transform);
			}

			printf("ENDMDL\n\n");
		}
	}

	JessQuery_free(Q);
	Molecule_free(M);
}

static Jess *init(const char *filename)
{
	FILE *file;
	FILE *temp;
	char buf[0x100];
	const char *s;
	Template *T;
	Jess *J;

	int line;
	int err;
	int k;

	if(!(file=fopen(filename,"r")))
	{
		perror(filename);
		exit(1);
	}

	J=Jess_create();
	line=0;
	while(fgets(buf,0x100,file))
	{
		line++;

		// Strip out blank lines and leading/trailing
		// spaces from the line...

		for(s=buf; isspace(*s); s++);
		for(k=strlen(s); k>0 && isspace(s[k-1]); k--);
		buf[k]=0;
		if(strlen(s)==0) continue;

		// Open the template file if possible.

		if(!(temp=fopen(s,"r")))
		{
			err=errno;
			fprintf(
				stderr,
				"%s, line %i: %s: %s\n",
				filename,
				line,
				s,
				strerror(err)
				);

			fclose(file);
			exit(1);
		}

		// Create a template from it if possible.

		if(!(T=TessTemplate_create(temp,s)))
		{
			fprintf(
				stderr,
				"%s, line %i: %s: error parsing template\n",
				filename,
				line,
				s
				);

			fclose(temp);
			fclose(file);
			exit(1);
		}
		fclose(temp);

		// Add it to the list of templates...

		Jess_addTemplate(J,T);
	}

	fclose(file);
	return J;
}

static void help(void)
{
	fprintf(
		stderr,
		"Jess version 0.4(gamma)\n"
		"Copyright (c) Jonathan Barker, 2002\n"
		"Command line syntax:\n\n"
		"   jess <T> <S> <r> <d> [F]\n\n"
		"where\n\n"
		"   <T> is the name of the template list file\n"
		"   <S> is a file containing a list of PDB filenames (use - for stdin)\n"
		"   <r> is the RMSD threshold\n"
		"   <d> is the distance cutoff\n"
		"   [F] are optional flags as a string with no spaces:\n"
		"         f: see PDB filenames in progress on stderr\n"
		"         n: do not transform coordinates of hit into\n"
		"            the template coordinate frame\n\n"
		"Contact jbarker@ebi.ac.uk for licensing\n"
		);

	exit(1);
}
// ==================================================================
// Entry point
// ==================================================================
// Arguments:
//	1				A file containing template filenames
//	2				A file containing PDB filenames
//	3				RMSD threshold (default 2)
//	4				Distance threshold (default 1)
// ==================================================================

int main(int argc, char **argv)
{
	FILE *file;
	char buf[0x100];
	const char *s;
	double tRmsd;
	double tDistance;
	Jess *J;
	int line,k;
	int count;
	//Riziotis edit
	int no_transform=0;

	if(argc<5 || argc>6) help();

	// Get optional flags

	if(argc==6)
	{
		for(s=argv[5]; *s; s++)
		{
			if(*s=='f') feedbackQ=1;
			//Riziotis edit
			else if(*s=='n') no_transform=1;
			else help();
		}
	}

	J=init(argv[1]);
	tRmsd=atof(argv[3]);
	tDistance=atof(argv[4]);

	if(strcmp(argv[2],"-")==0)
	{
		file=stdin;
	}
	else if(!(file=fopen(argv[2],"r")))
	{
		perror(argv[2]);
		exit(1);
	}

	line=0;
	while(fgets(buf,0x100,file))
	{
		line++;

		// Strip out blank lines and leading/trailing
		// spaces from the line...

		for(s=buf; isspace(*s); s++);
		for(k=strlen(s); k>0 && isspace(s[k-1]); k--);
		buf[k]=0;
		if(strlen(s)==0) continue;

		if(feedbackQ) fprintf(stderr,"%s\n",s);
		search(buf,J,tRmsd,tDistance,no_transform);
	}

	fclose(file);

	return 0;
}

// ==================================================================


