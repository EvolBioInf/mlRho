/***** formatPro.c ********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Oct 15 20:36:50 2012
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include "fp_interface.h"
#include "eprintf.h"
#include "fp_profileTree.h"

void scanFile(FILE *fp, Args *args);

int main(int argc, char *argv[]){
  int i, fd;
  char *version, *fileName;
  Args *args;
  FILE *fp;
  Profile *profiles;
  Node *root;
  ContigDescr *cd;

  version = "0.5";
  setprogname2("formatPro");
  args = getArgs(argc, argv);
  fileName = (char *)emalloc(256*sizeof(char));
  if(args->v)
    printSplash(version);
  if(args->h || args->e)
    printUsage(version);
  fileName = strcpy(fileName,args->n);
  fileName = strcat(fileName,".pos");
  fp = fopen(fileName,"wb");
  fwrite("pos",sizeof(char),3,fp);
  root = NULL;
  profiles = NULL;
  cd = (ContigDescr *)emalloc(sizeof(ContigDescr));
  cd->n = 0;
  cd->len = NULL;
  if(args->numInputFiles == 0){
    fd = 0;
    root = writePositions(fd, fp, root, cd, args);
  }else{
    for(i=0;i<args->numInputFiles;i++){
      fd = open(args->inputFiles[i], O_RDONLY, 0);
      root = writePositions(fd, fp, root, cd, args);
      close(fd);
    }
  }
  fclose(fp);
  profiles = getProfileArray();
  printf("#Positions written to %s\n", fileName);
  /* write profiles to file */
  writeProfiles(profiles, args);
  /* write contig lengths to file */
  cd = getContigDescr();
  writeContigs(cd, args);
  freeProfileArray(profiles);
  freeProfileTree(root);
  free(args);
  free(progname());
  free(fileName);
  return 0;
}

