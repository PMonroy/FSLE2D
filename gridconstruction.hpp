#ifndef GCONST
#define GCONST

#include "vectorXY.hpp"
#include "vectorIJ.hpp"
 
int MakeRegularGrid(vector<vectorXY> *points, vectorIJ *dimension, vectorXY domainmin, vectorXY intergrid, vectorXY domainmax);
vector<int> GetNeighbors(vectorIJ dimension);

#endif
