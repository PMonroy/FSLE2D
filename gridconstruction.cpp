#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#include "vectorXY.h"

int gridfsle2d(vector<vectorXY> *itracer, int *ni, int *nj, vectorXY domainmin, vectorXY intergrid, vectorXY domainmax)
{
  double x,y;
  int i,j;
  int ntracers;

  ntracers = (int)((domainmax.y-domainmin.y)/intergrid.y)+1;
  ntracers *= (int)((domainmax.x-domainmin.x)/intergrid.x)+1;

  if(ntracers<=0)
    return 1;
  
  (*itracer).reserve(ntracers);

  for(y=domainmin.y,j=0; y<=domainmax.y; y+=intergrid.x,j++)
    {
      for(x=domainmin.x,i=0; x<=domainmax.x; x+=intergrid.x,i++)
	{
	  (*itracer).push_back(vectorXY(x,y));
	}
    }

  *ni=i;
  *nj=j;

  return 0;
}

vector<int> neighbors(int ni, int nj)
{

  int i;
  int j;
  int ntracers=ni*nj;

  int q;
  int qmin[4],qincr[4],qmax[4];
  int qneighbor;
  int dir;

  vector<int> neighbor;

  neighbor.reserve(4*ntracers);

  qincr[0]=1;
  qincr[1]=ni;
  qincr[2]=-1;
  qincr[3]=-ni;

  qmin[1]=qmin[3]=0;
  qmax[1]=qmax[3]=ntracers;
 
  for(j=0; j<nj; j++)
    {
      qmin[0]=qmin[2]=j*ni;
      qmax[0]=qmax[2]=(j+1)*ni;
      for(i=0; i<ni; i++)
	{
	  q=i+j*ni;
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
