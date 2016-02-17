#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

#include "readparameters.h"
#include "gridconstruction.h" 
#include "constants.h"
#include "vflow.h" 
#include "integration.h"



int (*velocity)(double ,vectorXY , vectorXY* );

int main(int argc, char **argv)
{

  /**********************************************
   * READ COMAND LINE PARAMETERS
   **********************************************/
  string fnameparams;
  if(GetcmdlineParameters(argc, argv, &fnameparams))
    return 1;

  /* VERBOSE */
  if(verbose)
    {
      cout << "READ COMMAND LINE ARGUMENTS:" <<endl;
      cout << " Parameters file: " << fnameparams <<endl;
      cout << endl;
    }


  /**********************************************
   * READ PARAMETERS FROM FILE
   **********************************************/
  if(GetfileParameters(fnameparams))
    return 1;

  /* VERBOSE */
  if(verbose == 1)
    {
      cout << "PARAMETERS FROM FILE: "<< fnameparams <<endl; 
      cout << " vfield = "<<vfield<<endl ;
      cout << " domainmin = "<< domainmin<<endl;
      cout << " intergrid = " <<intergrid<<endl;
      cout << " domainmax = "<< domainmax<<endl;
      cout << " seeddate = "<< seeddate.tm_mday<<"-"<<seeddate.tm_mon+1<<"-"<<seeddate.tm_year<<endl ;
      cout << " intstep = "<<intstep<<endl ;
      cout << " tau = "<<tau<<endl ;
      cout << " deltamax = "<<deltamax<<endl ;
      cout << endl;
    }

  /**********************************************
   * GRID CONSTRUCTION
   **********************************************/
  vector<vectorXY> grid;
  int ni,nj;
  if(gridfsle2d(&grid, &ni, &nj, domainmin, intergrid, domainmax))
    {
      cout << "Error in contruction grid" << endl;
      return 1;
    }
  if(verbose == 1)
    {
      cout << "GRID CONSTRUCTION: "<<endl; 
      cout << " num. nodes = "<< grid.size()<<endl;
      cout << " ni = "<<ni<<endl;
      cout << " nj = "<<nj<<endl;
      cout << endl;
    }

  vector<int> neighbor;
  neighbor.reserve(4*grid.size());
  neighbor = neighbors(ni, nj);



  /**********************************************
   * INITIAL relative distances
   **********************************************/

  vectorXY delta,scalefactor;
  vector<double> ilength;
  unsigned int p;

  ilength.reserve(4*grid.size());

  for(unsigned int q=0; q<grid.size(); q++)
    {
      p=4*q;
      for(int dir=0; dir<4; dir++)
	{
	  if(neighbor[p+dir]>=0)
	    {
	      delta=grid[neighbor[p+dir]]-grid[q];
	      
	      delta.x=rads*delta.x;
	      delta.y=rads*delta.y;
	      
	      scalefactor.x=rearth*cos(rads*grid[q].y); 
	      scalefactor.y=rearth; 

	      delta=delta*scalefactor;
	      delta*=delta;

	      ilength.push_back(sqrt(delta.x+delta.y));
	    }
	  else
	    ilength.push_back(-1.0);
	}
    }

  /**********************************************
   * INITIALIZE VARIBLES
   **********************************************/

  vector<double> exit_time;
  vector<double> response;
  vector<int> qcore;
  vector<int> qflag;

  exit_time.reserve(grid.size());
  response.reserve(grid.size());
  qcore.reserve(grid.size());
  qflag.reserve(grid.size());


  for(unsigned int q=0; q<grid.size(); q++)
    {
      p=4*q;
      if(neighbor[p]!=-1 && neighbor[p+1]!=-1 && neighbor[p+2]!=-1 && neighbor[p+3]!=-1)
	{
	  qflag.push_back(1);// Initially we can compute fsle in this points (4=num. neighbours)
	  qcore.push_back(q);
	}
      else
	{
	  qflag.push_back(0); // Border points: we can not compute fsle on it (4>num. neighbours)
	}

      response.push_back(1.0);
      exit_time.push_back(1.0);
    }

  /**********************************************
   * TIME LOOP
   **********************************************/

  //Setup time parameters
  int ntime = abs(tau);
  struct tm *inidate = {0};
  double tstart;
  double tend;
  double h;
  int ascnd;

  ascnd = tau > 0;
  time_t seedtime = mktime(&seeddate); // convert date to time in seconds 
  if(ascnd)
    {
      time_t initime = seedtime;
      inidate = gmtime(&initime);
      tstart = 0.0;
      tend = (double) ntime;
      h=intstep;
    }
  else
    {
      time_t initime = seedtime - ntime*secondsday;
      inidate = gmtime(&initime); 
      tend = 0.0;
      tstart = (double) ntime;
      h=-1.0*intstep;
    }


  // Load velocity grid
  cout << "Loading grid:";
  if(loadvgrid(seeddate,vfield)!=0)
    {
      cout << " Error in reading grid from nc" << endl;
      return 1;
    }
  cout << " (Complete)"<<endl;

  // Load velocity field
  cout << "Loading velocities:";
  if(loadvflow(*inidate, ntime+2, vfield)!=0)
    {
      cout << "Error in reading velocities"<< endl;
      return 1;
    }
  cout << "(Complete)"<<endl;

  double t;
  int count;
  
  vector<vectorXY> tracer;
  tracer = grid;
  
  vector<double> length;
  length = ilength;
  
  double lengthmax;
  int dirmax;
  
  int q0, q1, q2, q3; 
  int qdir;

  cout << "Calculation of response and exit_time:" << endl;

  for(t=tstart,count=0; ascnd==1?(t<tend):(t>=tend); t+=h,count++)
     {

       cout << "count=" << count << "(" <<(tend-tstart)/h <<")"<<endl;

       /* Compute the position of tracer in time t=t+h*/
       for (unsigned int q = 0; q<tracer.size() ; q++)
	 {
	   //index of four q-neighboirs
	   q0 = neighbor[4*q];
	   q1 = neighbor[4*q+1];
	   q2 = neighbor[4*q+2];
	   q3 = neighbor[4*q+3];
	   
	   if(qflag[q]==-1)
	       continue;// skip the rest of the loop if q is a non-integrable point 
	     
	   if((qflag[q]=1) ||
	      (q0!=-1 && qflag[q0]=1) ||
	      (q1!=-1 && qflag[q1]=1) ||
	      (q2!=-1 && qflag[q2]=1) ||
	      (q3!=-1 && qflag[q3]=1))
	     {
	       if(RK4(t, h, &tracer[q], getvflow, vfield)==1)
		 {
		   qflag[q]=-1; // q is a non-integrable grid point
		   if(q0!=-1 && qflag[q0]!=-1) qflag[q0]=0;// We keep non-integrable neighbour grid points
		   if(q1!=-1 && qflag[q1]!=-1) qflag[q1]=0;
		   if(q2!=-1 && qflag[q2]!=-1) qflag[q2]=0;
		   if(q3!=-1 && qflag[q3]!=-1) qflag[q3]=0;
		 }
	     }
	 }
       /* Compute the relative distances */
       for(unsigned int q=0; q<tracer.size(); q++)
	 {
	   p=4*q;
	   for(int dir=0; dir<2; dir++)
	     {
	       qdir = neighbor[p+dir];
	       if((qflag[q]>0 && qdir>=0) || (qdir>=0 && qflag[qdir]>0))
		 {
		   delta=tracer[qdir]-tracer[q];
		   
		   delta.x=rads*delta.x;
		   delta.y=rads*delta.y;
		   
		   scalefactor.x=rearth*cos(rads*grid[q].y); 
		   scalefactor.y=rearth;
		   
		   delta=delta*scalefactor;
		   delta*=delta;
		   
		   length[p+dir]=sqrt(delta.x+delta.y);
		   length[4*qdir+2]=length[p+dir];
		 }
	     }
	 }

       /* Compute the max length*/
       for(unsigned int q=0; q<tracer.size(); q++)
	 {
	   if(qflag[q]>0)
	     {
	       p=4*q;
	       lengthmax=length[p];
	       for(int dir=0; dir<2; dir++)
		 {
		   if(length[p+dir]>lengthmax)
		     {
		       lengthmax=length[p+1];
		       dirmax=dir;
		     }
		   if(length[p+dir+2]>lengthmax)
		     {
		       lengthmax=length[p+dir+2];
		       dirmax=dir+2;
		     }
		 }

	       if(lengthmax>deltamax)
		 {
		   qflag[q]=0;
		   exit_time[q]=(count+1.0)*h;
		   response[q]=lengthmax/ilength[p+dirmax];
		 }
	     }
	 }
     }

  cout << "(Complete)" << endl;

  /****************
   * FREE MEMORY
   ****************/

  //freevgrid();
  //freevflow(ntime+2);

  /**********************************************************
   * COMPUTE FINITE SIZE LYAPUNOV EXPONENT
   **********************************************************/

  vector<double> fsle;
  fsle.reserve(grid.size());

  for(unsigned int q=0; q<qcore.size(); q++)
    {
      fsle.push_back(log(response[qcore[q]])/exit_time[qcore[q]]);
    }

 /**********************************************************
   * WRITING RESULT IN VTK FILE
   **********************************************************/
  ofstream ofile("fsle_debug.vtk");

  ofile<<"# vtk DataFile Version 3.0"<<endl;
  ofile<<"Complete vector field of ROMS Benguela"<<endl; 
  ofile<<"ASCII"<<endl;
  ofile<<"DATASET STRUCTURED_GRID"<<endl;
  ofile<<"DIMENSIONS "<< ni-2 <<" "<< nj-2 <<" "<< 1 <<endl;
  ofile<<"POINTS "<< (ni-2)*(nj-2) <<" float"<<endl;
  for(unsigned int q=0; q<qcore.size(); q++) 
    {
      ofile<<grid[qcore[q]].x<<" "<<grid[qcore[q]].y<<" "<< 0.0 <<endl;
    }
  ofile<<endl;
  ofile<<"POINT_DATA "<<(ni-2)*(nj-2)<<endl;
  ofile<<"SCALARS fsle float"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<fsle.size(); q++) 
    {
      ofile<<fsle[q]<<endl;
    }
  ofile<<"SCALARS qflag int"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<qcore.size(); q++) 
    {
      ofile<<qflag[qcore[q]]<<endl;
    }

  return 0;
}

