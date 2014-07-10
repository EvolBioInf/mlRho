/***** inspectPro.c *******************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Fri Oct 19 09:41:32 2012
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "ip_interface.h"
#include "ip_profileTree.h"
#include "eprintf.h"
#include "inspectPro.h"

void writeProfiles(Profile *pro, int n){
  int i;

  printf("#ID\tCount\tA\tC\tG\tT\n");
  for(i=0;i<n;i++)
    printf("%d\t%d\t%d\t%d\t%d\t%d\n",i,pro[i].n,
	   pro[i].profile[0],
	   pro[i].profile[1],
	   pro[i].profile[2],
	   pro[i].profile[3]);
}

void writeLengths(int *len, int n){
  int i;
  
  printf("#ID\tLength\n");
  for(i=0;i<n;i++)
    printf("%d\t%d\n",i,len[i]);
}

void writePositions(FILE *fp){
  int numRead;
  Position *pos;
  
  printf("#Pos\tPro\n");
  pos = (Position *)emalloc(sizeof(Position));
  numRead = fread(pos,sizeof(Position),1,fp);
  assert(numRead == 1);
  while(!feof(fp)){
    printf("%d\t%d\n",pos->pos,pos->pro);
    numRead = fread(pos,sizeof(Position),1,fp);
  }
  free(pos);
}

void writeLikelihoods(FILE *fp){
  int i;
  Result *r;
  double n, numRead, *lOne, *lTwo;

  r = (Result *)emalloc(sizeof(Result));
  printf("#id\tl_1\tl_2\n");
  numRead = fread(r,sizeof(Result),1,fp);
  assert(numRead == 1);
  numRead = fread(&n,sizeof(double),1,fp);
  lOne = (double *)emalloc(n*sizeof(double));
  lTwo = (double *)emalloc(n*sizeof(double));
  numRead = fread(lOne,sizeof(double),n,fp);
  assert(numRead == n);
  numRead = fread(lTwo,sizeof(double),n,fp);
  assert(numRead == n);
  for(i=0;i<n;i++)
    printf("%d\t%.3e\t%.3e\n",i,lOne[i],lTwo[i]);
  free(r);
  free(lOne);
  free(lTwo);
}

void scanFile(FILE *fp, Args *args){
  char *tag;
  int numRead, *lengths;
  Profile *profiles;
  int n;

  tag = (char *)emalloc(3*sizeof(char));
  numRead = fread(tag,sizeof(char),3,fp);
  assert(numRead == 3);
  while(!feof(fp)){
    if(strcmp(tag,"sum")==0){
      numRead = fread(&n,sizeof(int),1,fp);
      assert(numRead == 1);
      profiles = (Profile *)emalloc(n*sizeof(Profile));
      numRead = fread(profiles,sizeof(Profile),n,fp);
      assert(numRead == n);
      writeProfiles(profiles, n);
      free(profiles);
    }else if(strcmp(tag,"con")==0){
      numRead = fread(&n,sizeof(int),1,fp);
      assert(numRead == 1);
      lengths = (int *)emalloc(n*sizeof(int));
      numRead = fread(lengths,sizeof(int),n,fp);
      assert(numRead == n);
      writeLengths(lengths, n);
      free(lengths);
    }else if(strcmp(tag,"pos")==0){
      writePositions(fp);
    }else if(strcmp(tag,"lik")==0){
      writeLikelihoods(fp);
    }else
      assert(0);
    numRead = fread(tag,sizeof(char),3,fp);
  }
  free(tag);
}

int main(int argc, char *argv[]){
  int i;
  char *version;
  Args *args;
  FILE *fp;

  version = "0.3";
  setprogname2("inspectPro");
  args = getArgs(argc, argv);
  if(args->v)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);
  if(args->numInputFiles == 0){
    fp = stdin;
    scanFile(fp, args);
  }else{
    for(i=0;i<args->numInputFiles;i++){
      fp = efopen(args->inputFiles[i],"rb");
      scanFile(fp, args);
      fclose(fp);
    }
  }
  free(args);
  free(progname());
  return 0;
}


