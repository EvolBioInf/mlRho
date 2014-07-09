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

double numPos;
Node **resetProfilePairs(Node **profilePairs, int numProfiles);
Node **countPairs(Node **pairs, int numProfiles, char *inclPro, ContigDescr *contigDescr, FILE *fp, int dist, int memory);

Node **getProfilePairs(int numProfiles, char *inclPro, ContigDescr *contigDescr, FILE *fp, Args *args, int d){
	static Node **profilePairs = NULL;
	int i;

	numPos = 0;
	profilePairs = resetProfilePairs(profilePairs,numProfiles);
	if(args->L){
		for(i=0;i<args->S;i++)
			profilePairs = countPairs(profilePairs, numProfiles, inclPro, contigDescr, fp, d+i, args->o);
	}else
		profilePairs = countPairs(profilePairs, numProfiles, inclPro, contigDescr, fp, d, args->o);

	return profilePairs;
}

Node **resetProfilePairs(Node **profilePairs, int numProfiles){
	int i;

	if(profilePairs == NULL){
		profilePairs = (Node **)emalloc(numProfiles * sizeof(Node *));
		for(i=0;i<numProfiles;i++)
			profilePairs[i] = NULL;
	}else{
		for(i=0;i<numProfiles;i++){
			freeTree(profilePairs[i]);
			profilePairs[i] = NULL;
		}
	}
	return profilePairs;
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

Node **countPairs(Node **pairs, int numProfiles, char *inclPro, ContigDescr *contigDescr, FILE *fp, int dist, int memory){
	int a, b, i, l, r, tmp, numRead, bound;
	Position *pb;

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
					numPos++;
				}
			}
			r++;
			l++;
		}
	}
	return pairs;
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

double getNumPos(){
	return numPos;
}

