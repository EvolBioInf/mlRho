/***** ld.c ***************************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Thu Oct 25 08:50:39 2012
 **************************************************/
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "eprintf.h"
#include "interface.h"
#include "ld.h"

ContigDescr *thisContigDescr = NULL;

void setContigDescr(ContigDescr *contigDescr){
	thisContigDescr = contigDescr;
}


FILE *iniLdAna(Args *args){
	char *fileName;
	char *tag;
	FILE *fp;
	ContigDescr *cp;
	int i, max, n;

	int check;
	size_t numRead;

	/* get contig lengths from file */
	tag = (char *)emalloc(4*sizeof(char));
	
	check = asprintf( &fileName, "%s.con", args->n);
	assert( check != -1);

	fp = efopen(fileName,"rb");
	numRead = fread(tag,sizeof(char),3,fp);
	assert(numRead == 3);
	tag[3] = '\0';

	assert( strcmp(tag,"con") == 0 );

	cp = (ContigDescr *)emalloc(sizeof(ContigDescr));
	numRead = fread(&n,sizeof(int),1,fp);
	assert(numRead == 1);
	cp->pos = NULL;
	cp->n = n;
	cp->len = (int *)emalloc(n*sizeof(int));
	for(i=0;i<cp->n;i++){
		numRead = fread(&cp->len[i],sizeof(int),1,fp);
		assert(numRead == 1);
	}
	fclose(fp);
	max = 0;
	for(i=0;i<n;i++){
		if(cp->len[i] > max)
			max = cp->len[i];
	}
	cp->posBuf = (Position *)emalloc(max*sizeof(Position));
	setContigDescr(cp);

	free(fileName);
	/* get file pointer for position file */
	check = asprintf( &fileName, "%s.pos", args->n );
	assert( check != -1);

	fp = efopen(fileName,"rb");
	numRead = fread(tag,sizeof(char),3,fp);
	assert( strcmp(tag, "pos") == 0 );

	free(fileName);
	free(tag);
	return fp;
}

ContigDescr *getContigDescr(){
	return thisContigDescr;
}
