#include <iostream>
#include <iomanip>
#include <netcdfcpp.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>

using namespace std;

#include "vectorXY.hpp"
#include "vectorIJ.hpp"
#include "constants.hpp"
#include "rparameters.hpp"

static const int NC_ERR = 2;

//GLOBAL VARIABLES
unsigned int nvlon, nvlat;
unsigned int nv;
vectorIJ indexmin,indexmax;

vector <double> vlon; 
vector <double> vlat;
vectorXY ***vflow;

struct vflowParameters {

  const string dir;
  const int flat;

  const string slatdim;
  const string slondim;
  
  const string slatvar;
  const string slonvar;

  
  const double latstep;//must be in degrees
  const double lonstep;//must be in degrees
  
  const string suvar;
  const string svvar;

  const double uscalefactor;// factor to get velocity in meters per day
  const double vscalefactor;// factor to get velocity in meters per day

  const double ufillvalue;// velocity meters per day
  const double vfillvalue;// velocity meters per day


  // The constructor:
  vflowParameters(const string & vflowParamsFileName)
    :dir(getStringParam(vflowParamsFileName, "dir"))
    ,flat(getIntParam(vflowParamsFileName, "flat"))

    ,slatdim(getStringParam(vflowParamsFileName, "slatdim"))
    ,slondim(getStringParam(vflowParamsFileName, "slondim"))

    ,slatvar(getStringParam(vflowParamsFileName, "slatvar"))
    ,slonvar(getStringParam(vflowParamsFileName, "slonvar"))

    ,latstep(getDoubleParam(vflowParamsFileName, "latstep"))
    ,lonstep(getDoubleParam(vflowParamsFileName, "lonstep"))

    ,suvar(getStringParam(vflowParamsFileName, "suvar"))
    ,svvar(getStringParam(vflowParamsFileName, "svvar"))

    ,uscalefactor(getDoubleParam(vflowParamsFileName, "uscalefactor"))
    ,vscalefactor(getDoubleParam(vflowParamsFileName, "vscalefactor"))

    ,ufillvalue(getDoubleParam(vflowParamsFileName, "ufillvalue"))
    ,vfillvalue(getDoubleParam(vflowParamsFileName, "vfillvalue"))
{}

};
const vflowParameters vflowparams("vflow.conf");

//FUNCTIONS PROTOTYPES
vectorIJ GetIndices(vectorXY point);

//INLINE FUNCTIONS 
inline int LocateIndex(double x, double start, double step) {
  return floor((x-start)/step);
}
inline int LocateIndex(double x, const vector <double> &xx) {
  /* Given a vector xx[0...n-1], and given a value x, returns a value i such that x is between xx[i]
   * and xx[i+1]. xx must be monotonic, either increasing or decreasing. -1 or n is returned
   * to indicate that x is out of range.
   */

  if (x < xx.front()) 
    return -1;
  else if(x >= xx.back()) 
    return xx.size();
  else { 		
    unsigned long jl,jm,ju;
    jl=0;
    ju=xx.size()+1;
    int ascnd = (xx.back()>xx.front());
    while((ju-jl)>1) {
      jm =(ju+jl)>>1;
      if((x>=xx.at(jm-1))==ascnd) jl=jm;
      else ju=jm;
    }
    return int(jl-1);
  }
}

//FUNCTIONS DECLARATIONS
int LoadVGrid(struct tm rdate, vectorXY domainmin, vectorXY domainmax, vectorXY  meanvel, double duration) {
 
  NcError err(NcError::verbose_nonfatal);
  
  // Open the first Netcdf file
  char ncfile[256];
  sprintf(ncfile, "%s%04d-%02d-%02d.nc",vflowparams.dir.c_str(), rdate.tm_year,rdate.tm_mon+1,rdate.tm_mday);
  
  NcFile dataFile(ncfile, NcFile::ReadOnly);  
  if(!dataFile.is_valid()){
    cout <<ncfile <<endl;
    cout << "Error opening file:"<< ncfile <<endl;
    return NC_ERR;
  }

 // Get grid variable dimensions:
  NcDim *nvlonDim;
  NcDim *nvlatDim;

  if(!(nvlonDim = dataFile.get_dim(vflowparams.slondim.c_str()))
     || !(nvlatDim = dataFile.get_dim(vflowparams.slatdim.c_str()))) {
    cout << "Error getting grid variable dimensions"<<endl;
    return NC_ERR;
  }

  NcVar *vlonVar;
  NcVar *vlatVar;
  
  nvlon = nvlonDim->size();
  nvlat = nvlatDim->size();

  vlon.resize(nvlon);
  vlat.resize(nvlat);

  if (!(vlonVar = dataFile.get_var(vflowparams.slonvar.c_str()))
      || !vlonVar->get(&vlon[0], nvlon)){
    cout << "Error reading longitude variable"<<endl;
    return NC_ERR;
  }

  if (!(vlatVar = dataFile.get_var(vflowparams.slatvar.c_str()))
      || !vlatVar->get(&vlat[0], nvlat)){
    cout << "Error reading latitude variable"<<endl;
    return NC_ERR;
  }
  
#ifdef DEBUG
    cout << "Ncfile grid:" << endl;
    cout << " Longitude = "<< vlon.front()   <<","<< vlon.back()  <<" ("<< 0  <<","<< vlon.size()   <<")"<<endl;
    cout << " Latitude = " << vlat.front()   <<","<< vlat.back()  <<" ("<< 0  <<","<< vlat.size()   <<")"<<endl;
#endif

  vectorXY vdomainmin,vdomainmax;
  double scalefactor=(secondsday*duration*degrees)/rearth;

  vdomainmin=domainmin-meanvel*vectorXY(scalefactor/cos(rads*domainmin.y),scalefactor);
  vdomainmax=domainmax+meanvel*vectorXY(scalefactor/cos(rads*domainmax.y),scalefactor);

  indexmin=GetIndices(vdomainmin);
  indexmax=GetIndices(vdomainmax);
  indexmax+=vectorIJ(1,1);

  vlon.erase(vlon.begin()+indexmax.i,vlon.end());
  vlon.erase(vlon.begin(),vlon.begin()+indexmin.i);

  vlat.erase(vlat.begin()+indexmax.j,vlat.end());
  vlat.erase(vlat.begin(),vlat.begin()+indexmin.j);

  nvlon = vlon.size();
  nvlat = vlat.size();

  nv = nvlon*nvlat;
  
#ifdef DEBUG
    cout << "Velocity grid:" << endl;
    cout << " Longitude = "<< vlon.front()   <<","<< vlon.back()  <<" ("<< 0  <<","<< vlon.size()   <<")"<<endl;
    cout << " Latitude = " << vlat.front()   <<","<< vlat.back()  <<" ("<< 0  <<","<< vlat.size()   <<")"<<endl;
#endif

  return 0;
}
int LoadVFlow(struct tm seeddate, int ntime) {

  struct tm *date = {0};
  time_t seedtime, time; // Date in seconds UTC
  char ncfile[256];
  NcVar *uVar, *vVar;
  
  int t;
  
  vector <long int> vflowx_buffer(nv);
  vector <long int> vflowy_buffer(nv);
  
  vflow = new vectorXY **[ntime];
  
  // Get the Velocity data from the nc file
  seedtime = mktime(&seeddate); // convert date to time in seconds    
  for(t=0; t<ntime; t++){
    time = seedtime + t*secondsday;                                  
    date = gmtime(&time);
    
    sprintf(ncfile, "%s%04d-%02d-%02d.nc",vflowparams.dir.c_str(),date->tm_year, date->tm_mon+1, date->tm_mday);
    NcFile dataFile(ncfile, NcFile::ReadOnly);
    
#ifdef DEBUG
    //Reading nc file 
    cout << "Reading nc file: " << ncfile << " ";
#endif
    
    if(vflowparams.flat==0) {
      if(!dataFile.is_valid() 
	 || !(uVar = dataFile.get_var(vflowparams.suvar.c_str()))
	 || !(vVar = dataFile.get_var(vflowparams.svvar.c_str()))
	 || !uVar->set_cur(0, 0, indexmin.j, indexmin.i)
	 || !vVar->set_cur(0, 0, indexmin.j, indexmin.i)
	 || !uVar->get(&vflowx_buffer[0], 1, 1, nvlat, nvlon)
	 || !vVar->get(&vflowy_buffer[0], 1, 1, nvlat, nvlon)){
	cout << "[Fail]" <<endl;
	return NC_ERR;
      }
    }
    else if(vflowparams.flat==1){
      if(!dataFile.is_valid() 
	 || !(uVar = dataFile.get_var(vflowparams.suvar.c_str()))
	 || !(vVar = dataFile.get_var(vflowparams.svvar.c_str()))
	 || !uVar->set_cur(0, indexmin.j, indexmin.i)
	 || !vVar->set_cur(0, indexmin.j, indexmin.i)
	 || !uVar->get(&vflowx_buffer[0], 1, nvlat, nvlon)
	 || !vVar->get(&vflowy_buffer[0], 1, nvlat, nvlon)){
	cout << "[Fail]" <<endl;
	return NC_ERR;
      }
    }
    else{
      cout << "value of vflowparams.flat unknown:"<< vflowparams.flat <<endl;
      return 1;
    }
      
    unsigned int fj,q;
    
    vflow[t] = new vectorXY *[nvlat];      
    for(unsigned int j=0; j<nvlat; j++){
      fj=j*nvlon;
      vflow[t][j] = new vectorXY [nvlon];
      for(unsigned int i=0; i<nvlon; i++){
	q=fj+i;
	vflow[t][j][i].x= double(vflowx_buffer[q])*vflowparams.uscalefactor;
	vflow[t][j][i].y= double(vflowy_buffer[q])*vflowparams.vscalefactor;
      }
    }


#ifdef DEBUG
    //Success file read
    cout << "[Complete]" << endl;
#endif
  }

  return 0;
}

void FreeVFlow(unsigned int ntime) {

  unsigned int t, j;

  for(t=0; t<ntime; t++) {
    for(j=0; j<nvlat;j++) {
      delete[] vflow[t][j];
    }
    delete[] vflow[t];
  }
}

vectorIJ GetIndices(vectorXY point) {

  vectorIJ index;

  index.i = LocateIndex(point.x, vlon[0], vflowparams.lonstep);
  index.j = LocateIndex(point.y, vlat[0], vflowparams.latstep);

  if(index.i<0)
    index.i=0;
  else if((unsigned int) index.i>=vlon.size())
    index.i=vlon.size()-1;

  if(index.j<0)
    index.j=0;
  else if((unsigned int) index.j>=vlat.size())
    index.j=vlat.size()-1;

  return index;
}

int GetVFlow(double t,vectorXY point, vectorXY *vint) {

  vectorXY vgrid[8];
  vectorXY vcomp[8];
  vectorIJ index;
  unsigned int time;

  double alpha, beta;
  unsigned int i,j,l;
  unsigned int deltai,deltaj,deltatime;
  unsigned int q;

  time = (unsigned long) t;

  index.i = LocateIndex(point.x, vlon[0], vflowparams.lonstep);
  index.j = LocateIndex(point.y, vlat[0], vflowparams.latstep);

  /* Vectors and points with one flat index*/
  q=0;  
  for(deltatime=0; deltatime<2; deltatime++){
    l = time + deltatime;
    for(deltaj=0; deltaj<2; deltaj++){
      j = index.j + deltaj;
      for(deltai=0; deltai<2; deltai++){
	i = index.i + deltai;
	
	try{
	  vgrid[q].x = vlon.at(i);
	  vgrid[q].y = vlat.at(j);
	  }
	  catch(...){
	    return 1;
	  }
	
	vcomp[q] = vflow[l][j][i];	
	/* COAST CHECKING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
	if(vcomp[q].x <= vflowparams.ufillvalue || vcomp[q].y <= vflowparams.vfillvalue){
	  return 1;
	}
	q++; 
      }
    }
  }
  
  /* Longitude Interpolation: */
  alpha = (vgrid[1].x - point.x)/(vgrid[1].x - vgrid[0].x);
  beta = 1 - alpha;
  for(q=0; q<8; q+=2){
    vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
    vgrid[q>>1] = vgrid[q]; 
  }
  /* Latitude Interpolation: */
  alpha = (vgrid[1].y - point.y)/(vgrid[1].y - vgrid[0].y);
  beta = 1.0 - alpha;
  for(q = 0; q < 4; q+=2){
    vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
    vgrid[q>>1] = vgrid[q];
  }
  /* Time Interpolation: */ 
  alpha = ((double) (time + 1)) - t;  
  beta = 1.0 - alpha;   
  vcomp[0] = alpha * vcomp[0] + beta * vcomp[1];
  /* Interpolated V*/
  *vint = vcomp[0];

  return 0;
}
int GetVPartialDeriv(double t, vectorXY point, vectorIJ dir, vectorXY *dvdr){

  unsigned int time;
  unsigned int p;
  vectorXY delta;
  vectorIJ index;

  vectorXY vgrid[8];
  vectorXY partialderiv[8];
  double alpha,beta;
  unsigned int i,j,l;
  unsigned int i0,j0;
  unsigned int i1,j1;
  unsigned int deltai,deltaj,deltatime;

  index.i = LocateIndex(point.x, vlon[0], vflowparams.lonstep);
  index.j = LocateIndex(point.y, vlat[0], vflowparams.latstep);

  time = (unsigned int) t;

  p=0;
  for(deltatime=0; deltatime<2; deltatime++){
    l = time+deltatime;
    
    for(deltaj=0; deltaj<2; deltaj++){
      j  =  index.j + deltaj;
      if(j==0 || j==(nvlat-1))
	return 1;

      j0 = j - dir.j;
      j1 = j + dir.j;
      
      delta.x=rearth*cos(rads*vlat[j])*vflowparams.lonstep; 
      delta.y=rearth*vflowparams.latstep;
      
      for(deltai=0; deltai<2; deltai++) {
	
	i  = index.i + deltai;
	if(i==0 || i==(nvlon-1))
	  return 1;

	i0 = i - dir.i;
	i1 = i + dir.i;
	
	try{
	  vgrid[p].x = vlon.at(i);
	  vgrid[p].y = vlat.at(j);
	  
	  if(vgrid[p].x==vlon.front() 
	     || vgrid[p].y==vlon.back()
	     || vgrid[p].x==vlat.front() 
	     || vgrid[p].y==vlat.back()) {
	    return 1;
	  }
	}
	catch(...){
	  return 1;
	}
	
	//COASTAL CHECKING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(   vflow[l][j][i].x <= vflowparams.ufillvalue 
	   || vflow[l][j][i].y <= vflowparams.vfillvalue
	   || vflow[l][j0][i0].x <= vflowparams.ufillvalue 
	   || vflow[l][j0][i0].y <= vflowparams.vfillvalue
	   || vflow[l][j1][i1].x <= vflowparams.ufillvalue 
	   || vflow[l][j1][i1].y <= vflowparams.vfillvalue){
	  return 1;
	}
	
	partialderiv[p] = (vflow[l][j1][i1] - vflow[l][j0][i0])/(2.0*delta);
	p++;
      }
    }
  }
  
  /* Longitude Interpolation: */
  alpha = (vgrid[1].x - point.x)/(vgrid[1].x - vgrid[0].x);
  beta = 1 - alpha;
  for(p=0; p<8; p+=2){
    partialderiv[p>>1] = alpha * partialderiv[p] + beta * partialderiv[p+1];
    vgrid[p>>1] = vgrid[p]; 
  }

  /* Latitude Interpolation: */
  alpha = (vgrid[1].y - point.y)/(vgrid[1].y - vgrid[0].y);
  beta = 1.0 - alpha;
  for(p=0; p<4; p+=2){
    partialderiv[p>>1] = alpha * partialderiv[p] + beta * partialderiv[p+1];
    vgrid[p>>1] = vgrid[p];
  }
  
  /* Time Interpolation: */

  alpha = ((double) (time + 1)) - t;
  beta = 1 - alpha;

  partialderiv[0] = alpha*partialderiv[0]+beta*partialderiv[1];

  *dvdr = partialderiv[0];
  return 0;
}

