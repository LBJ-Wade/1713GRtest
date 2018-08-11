/*******************************************************************************
*                                                                              *
*  calcGalPdot.cc                                                              *
*                                                                              *
*  Calculates Pdot/P caused by differentia Galactic acceleration               *
*  for a pulsar at (l,b,d)                                                     *
*  Uses potential PJM11_best.Tpot
*                                                                              *
*  Based on C++ code written by Walter Dehnen, 1995-96,                        *
*                      Paul McMillan, 2007-,                                   *
*  Lund Observatory, Lund University.                                          *
*  address: Box 43, SE-221 00 Lund, Sweden                                     *
*  e-mail:  paul@astro.lu.se                                                   *
*                                                                              *
*******************************************************************************/

#include <fstream>
#include <iostream>
#include "GalPot.h"

using std::cout;
using std::cin;



int main(int argc,char *argv[])  
{
    ifstream file;
    
    double DEG2RAD = 0.017453292519943295;
    
    double Rsun=8.2,zsun=0.014;
    double xsun=-Rsun,ysun=0.0;
    
    double l,b,d;
    double nx,ny,nz;
    double R,x,y,z;
    double P,dPdR,dPdz;
    double gRsun,gxsun,gysun,gzsun,gR,gx,gy,gz;
    double gdiff,dDoD,dPoP;
    char inputfile[100];
    //char Model[15];
    
    if (argc != 3) {
        //cerr << "\n USAGE: surfacedensity.exe <R[kpc]> <z[kpc]> \n\n";
        cerr << "\n USAGE: surfacedensity.exe <Model> <R[kpc]> \n\n";
        exit(1);
    }
    
    //file.open("/Users/wex/Science/psrsoft/GalPot-McMillan2016/pot/PJM11_best.Tpot");
    //file.open("/home/zhuww/work/timing/GRtest/Lazaridis/GalPotMcMillan2016/pot/PJM16_best.Tpot");

    //strcpy(Model, argv[1]);
    sprintf(inputfile, "/homes/zhuww/data/timing/GRtest/Lazaridis/GalPotMcMillan2016/pot/%s.Tpot", argv[1]);
    //file.open("/homes/zhuww/data/timing/GRtest/Lazaridis/GalPotMcMillan2016/pot/PJM16_best.Tpot");
    file.open(inputfile);
    
    GalaxyPotential Phi(file);
    file.close();

	R = atof(argv[2]);
	//z = atof(argv[2]);

    //cout <<  R << " " << z << "\n";
    
    cout << Phi.DisksSurfaceDensity(R) << "\n" ;
    //cout << Phi.Mass(R) << "\n" ;


}
