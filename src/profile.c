/***** profile.c **********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Oct 24 12:24:34 2012
 **************************************************/
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "eprintf.h"
#include "profile.h"

void setNumProfiles(int numProfiles);
void setProfiles(Profile *profiles);

Profile *thisProfiles;
int thisNumProfiles;

void readProfiles(char *baseName){
	char *fileName, tag;
	int i, numProfiles, numRead;
	Profile *profiles;
	FILE *fp;

	int check = asprintf( &fileName, "%s.sum", baseName);

	fp = efopen(fileName,"rb");

	// FIXME: The following loop looks like a bug to me. Not only does it read the
	// same character thrice; Also the resulting `tag` is never used.
	for(i=0;i<3;i++)
		numRead = fread(&tag,sizeof(char),1,fp);
	assert(numRead == 1);

	numRead = fread(&numProfiles,sizeof(int),1,fp);
	assert(numRead == 1);

	profiles = (Profile *) emalloc(numProfiles * sizeof(Profile));
	for(i = 0; i < numProfiles; i++){
		numRead = fread( &profiles[i], sizeof(Profile), 1, fp);
		assert(numRead == 1);
	}

	setProfiles(profiles);
	setNumProfiles(numProfiles);
	fclose(fp);
	free(fileName);
}

void setProfiles(Profile *profiles){
	thisProfiles = profiles;
}

Profile *getProfiles(){
	return thisProfiles;
}

void setNumProfiles(int numProfiles){
	thisNumProfiles = numProfiles;
}

int getNumProfiles(){
	return thisNumProfiles;
}

double getCoverage(){
	long i, j;
	double c, sum, numProf;

	sum = 0;
	numProf = 0;
	for(i=0;i<thisNumProfiles;i++){
		c = 0;
		for(j=0;j<4;j++)
			c += thisProfiles[i].profile[j];
		sum += thisProfiles[i].n * c;
		numProf += thisProfiles[i].n;
	}
	return sum/numProf;
}
