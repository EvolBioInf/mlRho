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
#include "stringUtil.h"
#include "interface.h"
#include "profileTree.h"
#include "ld.h"

ProfilePairs *newProfilePairs(int numProfiles);
Node **resetProfilePairs(Node **profilePairs, int numProfiles);
ProfilePairs *countPairs(ProfilePairs *pp, char *inclPro, ContigDescr *contigDescr, FILE *fp, int dist, int memory);

ProfilePairs *getProfilePairs(int numProfiles, char *inclPro, ContigDescr *contigDescr, FILE *fp, Args *args, int d){
  int i;
  ProfilePairs *pp = NULL;

  pp = newProfilePairs(numProfiles);
  pp->numPos = 0;
  if(args->L){
    for(i=0;i<args->S;i++)
      pp = countPairs(pp, inclPro, contigDescr, fp, d+i, args->o);
  }else
    pp = countPairs(pp, inclPro, contigDescr, fp, d, args->o);

  return pp;
}

ProfilePairs *newProfilePairs(int numProfiles){
  int i;
  ProfilePairs *pp;

  pp = (ProfilePairs *)emalloc(sizeof(ProfilePairs));
  pp->numPos = 0;
  pp->numProfiles = numProfiles;
  pp->pairs = (Node **)emalloc(numProfiles*sizeof(Node *));
  for(i=0;i<numProfiles;i++)
    pp->pairs[i] = NULL;

  return pp;
}  

void freeProfilePairs(ProfilePairs *pp){
  int i;
  for(i=0;i<pp->numProfiles;i++)
    freeTree(pp->pairs[i]);
  free(pp->pairs);
  free(pp);
}

/* readPositions: read position information into contig descriptor
 *   assuming that fp points to the position file
 */
void readPositions(FILE *fp, ContigDescr *cd){
  int i, numRead;

  fseek(fp,3,SEEK_SET);
  cd->pos = (Position **)emalloc(cd->n*sizeof(Position *));
  for(i=0;i<cd->n;i++){
    cd->pos[i] = (Position *)emalloc(cd->len[i]*sizeof(Position));
    numRead = fread(cd->pos[i],sizeof(Position),cd->len[i],fp);
    assert(numRead == cd->len[i]);
  }

}

ProfilePairs *countPairs(ProfilePairs *pp, char *inclPro, ContigDescr *contigDescr, FILE *fp, int dist, int memory){
  int a, b, i, l, r, tmp, numRead, bound;
  Position *pb;
  Node **pairs;
 
  pairs = pp->pairs;
  pb = NULL;
  if(!memory){
    fseek(fp,3,SEEK_SET);
    pb = contigDescr->posBuf;
  }
  for(i=0;i<contigDescr->n;i++){
    if(memory)
      pb = contigDescr->pos[i];
    else{
      numRead = fread(pb,sizeof(Position),contigDescr->len[i],fp);
      assert(numRead == contigDescr->len[i]);
    }
    l = 0;
    r = 0;
    bound = contigDescr->len[i]-1;
    while(r<contigDescr->len[i]){
      while(pb[r].pos - pb[l].pos < dist && r < bound)
	r++;
      while(pb[r].pos - pb[l].pos > dist && l<=r)
  	l++;
      if(pb[r].pos-pb[l].pos == dist){
	a = pb[l].pro;
	b = pb[r].pro;
	if(inclPro[a] && inclPro[b]){
	  if(a > b){ /* sort indexes to halve the number of index pairs */
	    tmp = b;
	    b = a;
	    a = tmp;
	  }
	  pairs[b] = addTree(pairs[b],a);
	  pp->numPos++;
	}
      }
      r++;
      l++;
    }
  }
  return pp;
}

/* addTree: add key to tree */
Node *addTree(Node *node, int key){

  if(node == NULL)  /* new key has arrived */
    node = newNode(key);
  else if(node->key == key)
    node->n++;      /* repeated key */
  else if(node->key < key) /* descend into left subtree */
    node->left = addTree(node->left, key);
  else              /* descend into right subtree */
    node->right = addTree(node->right, key);

  return node;
}

/*newNode: generate and initialize new node */
Node *newNode(int key){
  Node *node;

  node = (Node *)emalloc(sizeof(Node));
  node->key = key;
  node->n = 1;
  node->left = NULL;
  node->right = NULL;

  return node;
}

/* freeTree: traverse tree and free the memory for each node */
void freeTree(Node *n){
  if(n != NULL){
    freeTree(n->left);
    freeTree(n->right);
    free(n);
  }
}


