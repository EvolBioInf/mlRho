/***** profileTree.c ******************************
 * Description: Tree for counting allele profiles.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Feb 18 17:06:35 2009.
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include "eprintf.h"
#include "fp_interface.h"
#include "fp_profileTree.h"

void setContigDescr(ContigDescr *contigDescr);

void writeIndexFile(int fd, Args *args);

char testMode = 0;
double numPos;
ContigDescr *globalContigDescr;
int foundNodeIndex;
int numNode = 0;
Node **globalNodeArray;

Node *writePositions(int fd, FILE *indF, Node *root, ContigDescr *contigDescr, Args *args){
  char *buf, *line;
  int n, i, j, l, maxNumInd, numInd, numWritten;
  short headerOpen;
  Node *node;
  Position position, *posArray;
  char **profileStr;

  node = (Node *)emalloc(sizeof(Node));
  node->key = (char *)emalloc(256*sizeof(char));
  node->c2 = 0;
  node->n = 0;
  node->id = 0;
  for(i=0;i<4;i++){
    node->profile1[i] = 0;
    node->profile2[i] = 0;
  }
  node->left = NULL;
  node->right = NULL;
  maxNumInd = 100000;
  profileStr = (char **)emalloc(4*sizeof(char *));
  for(i=0;i<4;i++)
    profileStr[i] = (char *)emalloc(256*sizeof(char));
  posArray = (Position *)emalloc(maxNumInd*sizeof(Position));
  buf = (char *)emalloc(args->b*sizeof(char));
  line = (char *)emalloc(100*sizeof(char));
  headerOpen = 0;
  l = 0;
  numInd = 0;
  while((n = read(fd, buf, args->b)) > 0){
    for(i=0; i<n; i++){
      if(buf[i] == '>'){
	headerOpen = 1;
	contigDescr->n++;
	contigDescr->len = (int *)erealloc(contigDescr->len,contigDescr->n*sizeof(int));
	contigDescr->len[contigDescr->n-1] = 0;
      }else if(headerOpen && buf[i] == '\n'){  /* reach end of header */
	headerOpen = 0;
      }else if(!headerOpen){
	if(buf[i] != '\n')
	  line[l++] = buf[i];
	else{                                  /* line filled */
	  line[l] = '\0';
	  l = 0;
	  sscanf(line,"%d\t%s\t%s\t%s\t%s", &(position.pos), profileStr[0],
		 profileStr[1], profileStr[2], profileStr[3]);
	  node->key[0] = '\0';
	  node->c1 = 0;
	  for(j=0;j<3;j++){
	    node->profile1[j] = atoi(profileStr[j]);
	    strcat(node->key,profileStr[j]);
	    strcat(node->key," ");
	    node->c1 += node->profile1[j];
	  }
	  node->profile1[j] = atoi(profileStr[j]);
	  strcat(node->key,profileStr[j]);
	  node->c1 += node->profile1[j];
	  if(node->c1 >= args->c){
	    numPos++;
	    contigDescr->len[contigDescr->n-1]++;
	    root = addTree(root,node->key,1,node,NULL);
	    position.pro = foundNodeIndex;
	    posArray[numInd++] = position;
	    if(numInd == maxNumInd){
	      numWritten = fwrite(posArray,sizeof(Position),numInd,indF);
	      assert(numWritten == numInd);
	      numInd = 0;
	    }
	  }
	}
      }
    }
  }
  numWritten = fwrite(posArray,sizeof(Position),numInd,indF);
  setContigDescr(contigDescr);
  free(buf);
  free(line);
  free(node->key);
  free(node);
  free(posArray);
  for(i=0;i<4;i++)
    free(profileStr[i]);
  free(profileStr);

  return root;
}

void  writeProfiles(Profile *profiles, Args *args){
  char *fileName;
  FILE *profileF;
  /* int i; */

  if(profiles != NULL){
    fileName = (char *)emalloc(256*sizeof(char));
    fileName = strcpy(fileName,args->n);
    fileName = strcat(fileName,".sum");
    profileF = efopen(fileName,"wb");
    fwrite("sum",sizeof(char),3,profileF);
    fwrite(&numNode,sizeof(int),1,profileF);
    fwrite(profiles,sizeof(Profile),numNode,profileF);
    fclose(profileF);
    printf("#Profile summary written to %s\n", fileName);
  }else{
    printf("ERROR: something is wrong with the profile tree.\n");
    exit(-1);
  }
  free(fileName);
}

void writeContigs(ContigDescr *cd, Args *args){
  char *fileName;
  FILE *contF;

  fileName = (char *)emalloc(256*sizeof(char));
  fileName = strcpy(fileName,args->n);
  fileName = strcat(fileName,".con");
  contF = efopen(fileName,"wb");
  fwrite("con",sizeof(char),3,contF);
  fwrite(&cd->n,sizeof(int),1,contF);
  fwrite(cd->len,sizeof(int),cd->n,contF);
  printf("#Contig lengths written to %s\n",fileName);
  free(fileName);
  fclose(contF);
}

void setContigDescr(ContigDescr *contigDescr){
  globalContigDescr = contigDescr;
}

ContigDescr *getContigDescr(){
  return globalContigDescr;
}

/* addTree: add key to tree */
Node *addTree(Node *node, char *key, int count, Node *n1, Node *n2){
  int cond;

  if(node == NULL){  /* new key has arrived */
    node = newNode(key,count,n1,n2);
    foundNodeIndex = node->id;
  }else if((cond = strcmp(key, node->key)) == 0){
    node->n++;      /* repeated key */
    foundNodeIndex = node->id;
  }else if(cond < 0) /* descend into left subtree */
    node->left = addTree(node->left, key, count, n1, n2);
  else              /* descend into right subtree */
    node->right = addTree(node->right, key, count, n1, n2);

  return node;
}

/*newNode: generate and initialize new node */
Node *newNode(char *key, int count, Node *n1, Node *n2){
  Node *node;
  int i, l;
  static int id = 0;
  static Node **nodeArray = NULL;

  node = (Node *)emalloc(sizeof(Node));
  node->id = id++;
  l = strlen(key);
  node->key = emalloc(sizeof(char)*(l+1));
  strcpy(node->key,key);
  node->n = count;
  node->left = NULL;
  node->right = NULL;
  for(i=0;i<4;i++)
    node->profile1[i] = n1->profile1[i];
  node->c1 = n1->c1;
  if(n2){
    for(i=0;i<4;i++)
      node->profile2[i] = n2->profile1[i];
    node->c2 = n2->c1;
  }else{
    node->c2 = 0;
    for(i=0;i<4;i++)
      node->profile2[i] = 0;
  }
  numNode++;
  nodeArray = (Node **)erealloc(nodeArray,id*sizeof(Node *));
  nodeArray[id-1] = node;
  globalNodeArray = nodeArray;
  return node;
}

/* printTree: traverse tree and print profile & count for every node */
void printTree(FILE *fp, Node *node){

  if(node != NULL){
    printTree(fp, node->left);

    fprintf(fp,"%d|%s|%d %d %d %d",node->n,node->key,node->profile1[0],node->profile1[1],node->profile1[2],node->profile1[3]);
    if(node->profile2)
      fprintf(fp," %d %d %d %d",node->profile2[0],node->profile2[1],node->profile2[2],node->profile2[3]);
    fprintf(fp,"\n");

    printTree(fp, node->right);
  }
}

/* freeTree: traverse tree and free the memory for each node */
void freeTree(Node *n){
  if(n != NULL){
    freeTree(n->left);
    freeTree(n->right);
    free(n->key);
    free(n);
  }
}

void setTestMode(){
  testMode  = 1;
}

double getNumPos(){
  return numPos;
}

void freeProfileTree(Node *root){
  free(globalContigDescr->len);
  free(globalContigDescr);
  free(globalNodeArray);
  freeTree(root);
}

void freeProfileArray(Profile *pro){
  free(pro);
}

Profile *getProfileArray(){
  Profile *pa;
  int i, j;

  pa = (Profile *)emalloc(numNode*sizeof(Profile));
  for(i=0;i<numNode;i++){
    for(j=0;j<4;j++)
      pa[i].profile[j] = globalNodeArray[i]->profile1[j];
    pa[i].n = globalNodeArray[i]->n;
  }

  return pa;
}

