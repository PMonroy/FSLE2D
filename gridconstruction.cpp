#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

#include "vectorXY.hpp"
#include "vectorIJ.hpp"

int MakeRegularGrid(vector<vectorXY> *point, vectorIJ *dimension, vectorXY domainmin, vectorXY intergrid, vectorXY domainmax)
{
  double x,y;
  int i=0,j=0;
  int numpoints;

  numpoints = (int)((domainmax.y-domainmin.y)/intergrid.y)+1;
  numpoints *= (int)((domainmax.x-domainmin.x)/intergrid.x)+1;

  if(numpoints<1)
    return 1;
  
  (*point).reserve(numpoints);

  for(y=domainmin.y,j=0; y<=domainmax.y; y+=intergrid.x,j++){
    for(x=domainmin.x,i=0; x<=domainmax.x; x+=intergrid.x,i++){
	  (*point).push_back(vectorXY(x,y));
    }
  }

  //Dimension
  *dimension=vectorIJ(i,j);

  return 0;
}

vector<int> GetNeighbors(vectorIJ dimension) {

  int i;
  int j;
  int numpoints=dimension.i*dimension.j;
  
  int q;
  int qmin[4],qincr[4],qmax[4];
  int qneighbor;
  int dir;

  vector<int> neighbor;

  neighbor.reserve(4*numpoints);

  qincr[0]=1;
  qincr[1]=dimension.i;
  qincr[2]=-1;
  qincr[3]=-dimension.i;

  qmin[1]=qmin[3]=0;
  qmax[1]=qmax[3]=numpoints;
 
  for(j=0; j<dimension.j; j++)
    {
      qmin[0]=qmin[2]=j*dimension.i;
      qmax[0]=qmax[2]=(j+1)*dimension.i;
      for(i=0; i<dimension.i; i++)
	{
	  q=i+j*dimension.i;
	  for(dir=0; dir<4; dir++)
	    {
	      qneighbor=q+qincr[dir];
	      if(qneighbor>=qmin[dir] && qneighbor<qmax[dir])
		neighbor.push_back(qneighbor);
	      else
		neighbor.push_back(-1);
	    }
	}
    }

  return neighbor;
}
