#ifndef GCONST
#define GCONST

#include "vectorXY.h"
 
int gridfsle2d(vector<vectorXY> *itracer,int *ni, int *nj, vectorXY domainmin, vectorXY intergrid, vectorXY domainmax);

vector<int> neighbors(int ni, int nj);

#endif
