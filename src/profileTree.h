/***** profileTree.h ******************************
 * Description: Header file for profileTree.c.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Wed Feb 18 17:06:39 2009.
 * Distributed under the GNU General Public License
 **************************************************/
#ifndef PROFILETREE
#define PROFILETREE
#include <stdio.h>
#include "interface.h"
#include "ld.h"

typedef struct node{  /* the tree node: */
  int key;            /* sort key */
  int n;              /* number of occurrences */
  struct node *left;  /* left child */
  struct node *right; /* right child */
}Node;

typedef struct profilePairs{
  int numProfiles;
  int dist;
  long numPos;
  Node **pairs;
  double pi;
}ProfilePairs;

ProfilePairs *getProfilePairs(int numProfiles, char *inclPro, ContigDescr *contigDescr, FILE *fp, Args *args, int d);
void printTree(FILE *fp, Node *node);
void freeTree(Node *n);
void setTestMode();
double getNumPos();
int getNumNode();
double getNumPosBam();
Node *newNode(int key);
Node *addTree(Node *node, int key);
int coverage(char *key, int d);
void freeProfileTree();
Node *getSummarizedProfiles(Args *args);
FILE *iniLinkAna(Args *args);
int getNumNode();
void readPositions(FILE *fp, ContigDescr *cd);
void freeProfilePairs(ProfilePairs *pp);
#endif
