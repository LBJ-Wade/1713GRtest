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
    double gdiff_hor,dDoD_hor,dPoP_hor;
    double gdiff_ver,dDoD_ver,dPoP_ver;
    char inputfile[100];
    
    if (argc != 5) {
        cerr << "\n USAGE: calcGalPdot.exe <l[deg]> <b[deg]> <d[kpc]> <Model>\n\n";
        exit(1);
    }
    
    sprintf(inputfile, "/homes/zhuww/data/timing/GRtest/Lazaridis/GalPotMcMillan2016/pot/%s.Tpot", argv[4]);
    //file.open("/homes/zhuww/data/timing/GRtest/Lazaridis/GalPotMcMillan2016/pot/PJM16_best.Tpot");
    file.open(inputfile);
    //file.open("/Users/wex/Science/psrsoft/GalPot-McMillan2016/pot/PJM11_best.Tpot");
    //file.open("/home/zhuww/work/timing/GRtest/Lazaridis/GalPotMcMillan2016/pot/PJM16_best.Tpot");
    //file.open("/homes/zhuww/data/timing/GRtest/Lazaridis/GalPotMcMillan2016/pot/PJM16_best.Tpot");
    //file.open(inputfile);
    
    GalaxyPotential Phi(file);
    file.close();

	//l = std::stod(argv[1]);
	//b = std::stod(argv[2]);
	//d = std::stod(argv[3]);
	l = atof(argv[1]);
	b = atof(argv[2]);
	d = atof(argv[3]);
    
    /*
     * position of the pulsar in the Galaxy 
     */
    
    //unit vector (towards pulsar) 
	nx = cos(b*DEG2RAD) * cos(l*DEG2RAD);
	ny = cos(b*DEG2RAD) * sin(l*DEG2RAD);
	nz = sin(b*DEG2RAD);
	
	x = xsun + d*nx;
	y = ysun + d*ny;
	z = zsun + d*nz;
	
	R = sqrt(x*x + y*y);
	
	/*
	 * Galactic acceleration of the Sun
     */
     
    P = Phi(Rsun,zsun,dPdR,dPdz); //[kpc, kpc, kpc/Myr^2, kpc/Myr^2]
    
    gRsun = -dPdR * 3.09843657844953E-6; //[cm/s^2]
    gzsun = -dPdz * 3.09843657844953E-6; //[cm/s^2]
    
    gxsun = -gRsun;
    gysun =  0.0;
         
	/*
	 * Galactic acceleration of the Pulsar
     */

    P = Phi(R,z,dPdR,dPdz); //[kpc, kpc, kpc/Myr^2, kpc/Myr^2]
    
    gR = -dPdR * 3.09843657844953E-6; //[cm/s^2]
    gz = -dPdz * 3.09843657844953E-6; //[cm/s^2]
    
    gx = gR*x/R;
    gy = gR*y/R;
    
    /*
     * Differential acceleration
     */
     
     gdiff_hor = (gx - gxsun)*nx + (gy - gysun)*ny;
     gdiff_ver = (gz - gzsun)*nz;
     dDoD_hor  = -gdiff_hor / 2.99792458E10;
     dPoP_hor  = -dDoD_hor;
     dDoD_ver  = -gdiff_ver / 2.99792458E10;
     dPoP_ver  = -dDoD_ver;
     
    /*
     * Output
     */
    
    cout << dPoP_hor << " "  << dPoP_ver << "\n";
    //cout << (gx - gxsun) << (gy - gysun) << (gz - gzsun) << "\n";
}
