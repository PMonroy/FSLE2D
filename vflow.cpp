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

#include "vectorXY.h"
#include "constants.h"

static const int NC_ERR = 2;

struct  vectorIJ {
    int i;
    int j;
};
//EXTERN VARIABLES
extern int verbose;

//GLOBAL VARIABLES
int nvlon, nvlat;
vector <double> vlon; 
vector <double> vlat;
vector<vector<vector<vectorXY> > > vflow;

// FLOW MODEL PARAMETERS
enum vmodels { MYOCEAN,
	       AVISO,
	       NUMMODELS
};
struct velocitymodel {  
  string dir;
  string slatdim;
  string slondim;
  string slatvar;
  string slonvar;
  double latstep;//must be in degrees
  double lonstep;//must be in degrees
  string svvar;
  string suvar;
  double vscalefactor;// factor to get velocity in meters per day
  double uscalefactor;// factor to get velocity in meters per day
  double vfillvalue;// velocity meters per day
  double ufillvalue;// velocity meters per day
};
velocitymodel vmodel[NUMMODELS];

void loadparamsmodel(void )
{
  /**************************************
   * MYOCEAN PARAMETERS
   *************************************/
  vmodel[MYOCEAN].dir="/net/argos/data/peps/dflod/data/FSLEz/MODEL_FORECAST/netcdf/";
  vmodel[MYOCEAN].slatdim="latitude";
  vmodel[MYOCEAN].slondim="longitude";
  vmodel[MYOCEAN].slatvar="latitude";
  vmodel[MYOCEAN].slonvar="longitude";
  vmodel[MYOCEAN].latstep=1.0/12.0;//must be in degrees
  vmodel[MYOCEAN].lonstep=1.0/12.0;//must be in degrees
  vmodel[MYOCEAN].svvar="v";
  vmodel[MYOCEAN].suvar="u";
  
  vmodel[MYOCEAN].vscalefactor=0.000610370188951492;
  vmodel[MYOCEAN].vscalefactor*=secondsday;// factor to get velocity in meters per day
 
  vmodel[MYOCEAN].uscalefactor=0.000610370188951492;
  vmodel[MYOCEAN].uscalefactor*=secondsday;// factor to get velocity in meters per day

  vmodel[MYOCEAN].vfillvalue= -32767.0;
  vmodel[MYOCEAN].vfillvalue*=vmodel[MYOCEAN].vscalefactor;// velocity meters per day

  vmodel[MYOCEAN].ufillvalue= -32767.0;
  vmodel[MYOCEAN].ufillvalue*=vmodel[MYOCEAN].uscalefactor;// velocity meters per day

  /**************************************
   * AVISO PARAMETERS
   *************************************/
  vmodel[AVISO].dir="/net/argos/data/peps/dflod/data/FSLEz/AVISO/";
  vmodel[AVISO].slatdim="lat";
  vmodel[AVISO].slondim="lon";
  vmodel[AVISO].slatvar="lat";
  vmodel[AVISO].slonvar="lon";
  vmodel[AVISO].latstep=0.25;//must be in degrees
  vmodel[AVISO].lonstep=0.25;//must be in degrees
  vmodel[AVISO].svvar="v";
  vmodel[AVISO].suvar="u";

  vmodel[AVISO].vscalefactor= 0.0001;
  vmodel[AVISO].vscalefactor*=secondsday;// factor to get velocity in meters per day

  vmodel[AVISO].uscalefactor= 0.0001;
  vmodel[AVISO].uscalefactor*=secondsday;// factor to get velocity in meters per day

  vmodel[AVISO].vfillvalue= -2147483647.0;
  vmodel[AVISO].vfillvalue*=vmodel[AVISO].vscalefactor;// velocity meters per day
 
  vmodel[AVISO].ufillvalue= -2147483647.0;
  vmodel[AVISO].ufillvalue*=vmodel[AVISO].uscalefactor;// velocity meters per day
 
}

int loadvgrid(struct tm rdate, int vfield)
{
  
  char ncfile[256];
  NcError err(NcError::verbose_nonfatal);

  //Loading parameters of the model
  loadparamsmodel();

  // Open the first Netcdf file
  sprintf(ncfile, "%s%04d-%02d-%02d.nc",vmodel[vfield].dir.c_str(), rdate.tm_year,rdate.tm_mon+1,rdate.tm_mday);
  NcFile dataFile(ncfile, NcFile::ReadOnly);  

  // Check to see if the file was opened.
  if(!dataFile.is_valid())
    return NC_ERR;

  NcDim *nvlonDim;
  if (!(nvlonDim = dataFile.get_dim(vmodel[vfield].slondim.c_str())))
    return NC_ERR;
  nvlon = nvlonDim->size();

  NcDim *nvlatDim;
  if (!(nvlatDim = dataFile.get_dim(vmodel[vfield].slatdim.c_str())))
    return NC_ERR;
  nvlat = nvlatDim->size();

  NcVar *vlonVar;
  if (!(vlonVar = dataFile.get_var(vmodel[vfield].slonvar.c_str())))
    return NC_ERR;

  NcVar *vlatVar;
  if (!(vlatVar = dataFile.get_var(vmodel[vfield].slatdim.c_str())))
    return NC_ERR;

  vlon.resize(nvlon);
  vlat.resize(nvlat);

  // Get the lat/lon data from the nc file.
  if (!vlonVar->get(&vlon[0], nvlon))
    return NC_ERR;

  if (!vlatVar->get(&vlat[0], nvlat))
    return NC_ERR;

  return 0;
}

int loadvflow(struct tm seeddate, int ntime, int vfield)
{

  struct tm *date = {0};
  time_t seedtime, time; // Date in seconds UTC
  char ncfile[256];
  NcVar *uVar, *vVar;

  int i,j,k;
  int t;

  vector <long int> ubuffer;
  vector <long int> vbuffer;

  int nv = nvlat*nvlon;

  ubuffer.resize(nv);
  vbuffer.resize(nv);

  // Dynamically allocate Velocity field
  vflow.resize(ntime);
  for(k = 0; k < ntime; k++) 
    {
      vflow[k].resize(nvlat);
      for(j = 0; j < nvlat; j++) 
	{
	  vflow[k][j].resize(nvlon);
	}
    }

  // Get the Velocity data from the nc file
  seedtime = mktime(&seeddate); // convert date to time in seconds    
  for(t=0; t<ntime; t++)  
    {
      time = seedtime + t*secondsday;                                  
      date = gmtime(&time);

      sprintf(ncfile, "%s%04d-%02d-%02d.nc",vmodel[vfield].dir.c_str(),date->tm_year, date->tm_mon+1, date->tm_mday);
      NcFile dataFile(ncfile, NcFile::ReadOnly);

      if(verbose==1)
	cout << "Reading nc file: " << ncfile <<endl;

      // Check to see if the file was opened.
      if(!dataFile.is_valid())
	return NC_ERR;

      if (!(uVar = dataFile.get_var(vmodel[vfield].suvar.c_str())))
	return NC_ERR;

      if (!(vVar = dataFile.get_var(vmodel[vfield].svvar.c_str())))
	return NC_ERR;
      if(vfield == 0)
	{
	  if (!uVar->set_cur(0, 0, 0, 0))
	    return NC_ERR;
	  if (!vVar->set_cur(0, 0, 0, 0))
	    return NC_ERR;
	  
	  if (!uVar->get(&ubuffer[0], 1, 1, nvlat, nvlon))
	    return NC_ERR;
	  if (!vVar->get(&vbuffer[0], 1, 1, nvlat, nvlon))
	    return NC_ERR;
	}
      else
	{
	  if (!uVar->set_cur(0, 0, 0))
	    return NC_ERR;
	  if (!vVar->set_cur(0, 0, 0))
	    return NC_ERR;
	  
	  if (!uVar->get(&ubuffer[0], 1, nvlat, nvlon))
	    return NC_ERR;
	  if (!vVar->get(&vbuffer[0], 1, nvlat, nvlon))
	    return NC_ERR;
	}
      for(int q=0; q<nv; q++)
	{
	  j=(int) (q/nvlon);
	  i=(q-j*nvlon);
	  vflow[t][j][i]=vectorXY(((double) ubuffer[q])*vmodel[vfield].uscalefactor,
				  ((double) vbuffer[q])*vmodel[vfield].vscalefactor);
	}
    }

  /*  ofstream vfile("velocities.vtk");
  
  vfile<<"# vtk DataFile Version 3.0"<<endl;
  vfile<<"Complete data GLORY"<<endl; 
  vfile<<"ASCII"<<endl;
  vfile<<"DATASET STRUCTURED_GRID"<<endl;
  vfile<<"DIMENSIONS "<< nvlon <<" "<< nvlat <<" "<<1<<endl;
  vfile<<"POINTS "<< nvlon*nvlat <<" float"<<endl;
  
  for(int j=0;j<nvlat;j++) 
    {
      for(int i=0;i<nvlon;i++) 
	{
	  vfile << vlon[i] <<" "<< vlat[j]<<" "<< 0.0 <<endl;
	}
    }
  vfile<<endl;
  vfile<<"POINT_DATA "<< nvlon*nvlat <<endl;
  vfile<<"VECTORS velocity float"<<endl;
  for(int j=0;j<nvlat;j++) 
    {
      for(int i=0;i<nvlon;i++)
	{	
	  vfile<<vflow[ntime-1][j][i].x<<" ";
	  vfile<<vflow[ntime-1][j][i].y<<" ";
	  vfile<<0.0<<endl;
	}
    }
    vfile.close();*/
    
  return 0;
}

int GetIndices(vectorXY point, vectorIJ *index, int vfield)
{
  /* Locate index longitude*/
  index->i = floor((point.x-vlon[0])/vmodel[vfield].lonstep);

  if(index->i < 0 || index->i >= (nvlon-1))
      return 1;

  /* Locate index latitude*/
  index->j = floor((point.y-vlat[0])/vmodel[vfield].latstep);

  if(index->j < 0 || index->j >= (nvlat-1))
      return 1;

  return 0;
}

int getvflow(double t,vectorXY point, vectorXY *vint, int vfield)
{
  vectorXY vgrid[8];
  vectorXY vcomp[8];
  vectorIJ index;
  unsigned long time;

  double alpha, beta;

  int i,j,k;
  int deltatime,deltai,deltaj;
  unsigned int q;

  /* Calculates the vectors vb[15] and points ptmb[15] */
  time = (unsigned long) t;

  if(GetIndices(point, &index, vfield)==1)
      return 1;

  /* Vectors and points with only one int index*/
  q=0;  
  for(deltatime=0; deltatime<2; deltatime++)
    {
      for(deltai=0; deltai<2; deltai++)
	{
	  for(deltaj = 0; deltaj<2; deltaj++)
	    {
	      i = index.i + deltai;
	      j = index.j + deltaj;
	      k = time+deltatime;
	      
	      vgrid[q].x = vlon[i];
	      vgrid[q].y = vlat[j];
	      vcomp[q] = vflow[k][j][i];

	      /* COAST CHECKING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
	      if(vcomp[q].x == vmodel[vfield].ufillvalue || vcomp[q].y == vmodel[vfield].vfillvalue)
		  return 1;

	      q++; 
	    }
	}
    }
    


  /* Latitude Interpolation: */
  alpha = (vgrid[1].y - point.y)/(vgrid[1].y - vgrid[0].y);
  beta = 1.0 - alpha;
  for(q = 0; q < 8; q+=2)
    {
      vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
      vgrid[q>>1] = vgrid[q];
    }

  /* Longitude Interpolation */
  alpha = (vgrid[1].x - point.x)/(vgrid[1].x - vgrid[0].x);
  beta = 1.0 - alpha;
  for(q = 0; q < 4; q+=2)
    {
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
