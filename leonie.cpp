#include <stdio.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
//#include <unistd.h>
#include <math.h>
#include "dish.h"
#include "random.h"
#include "cell.h"
#include "parameter.h"
#include "sqr.h"


void PDE::Secrete(CellularPotts *cpm) {

  const double dt=par.dt;

  for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++) {
      // inside cells
      if (cpm->Sigma(x,y)) {

	sigma[0][x][y]+=par.secr_rate[0]*dt;

      } else {
      // outside cells
	sigma[0][x][y]-=par.decay_rate[0]*dt*sigma[0][x][y];

      }
    }
}


void PDE::InitializeAgeLayer(int l,double value,CellularPotts *cpm){
  for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++){
	if(sigma[l][x][y]>0){
	    sigma[l][x][y]=value;}
	else {sigma[l][x][y]=0.;}

    }
}

void CellularPotts::ChangeThetas(Dish *dish){
  int number_of_cells=dish->CountCells();
  double dt=0.1*par.tau;
  for (int s=1;s<number_of_cells+1;s++){
    dish->getCell(s).SetTheta(dish->getCell(s).theta + sqrt(2.0/par.tau)*generateGaussianNoise(0, sqrt(2*dt/par.tau)));
  }
}


void PDE::AgeLayer(int l,double value,CellularPotts *cpm, Dish *dish){
for (const auto& elem: cpm->alivePixels) {
    /* ... process elem ... */
    int x= (int)(elem/100000LL);
    int y= (int)(elem%100000LL);
    if(sigma[l][x][y]>0.)
      sigma[l][x][y]-=(value + 0.001*cpm->concentration[x][y]);

    // // do not allow young lattice sites outside the T cells
    // if (dish->getCell(cpm->Sigma(x,y)).getTau() ==0 )
    //   sigma[l][x][y]= 0.;
  }
}

void PDE::InitializeMILayer(int l,double value,CellularPotts *cpm){
  for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++){
 //cout << "Anything..." << endl;
    if (cpm->Sigma(x,y)==0 ){
	sigma[l][x][y]=0;
	//cout << cpm->Sigma(x,y) << " has now value " << sigma[l][x][y] << endl;
}
    else{
	//cout << cpm->Sigma(x,y);
	int rand_value=RandomNumber(par.max_matrix/4);
	sigma[l][x][y]=rand_value;
	//cout << " has now value " << sigma[l][x][y] << endl;
}

    }
}


void PDE::MILayer(int l,double value,CellularPotts *cpm, Dish *dish){
for (int x=0;x<sizex;x++)
    for (int y=0;y<sizey;y++){

    if(sigma[l][x][y]>=0.&&sigma[l][x][y]<par.max_matrix-1)
      sigma[l][x][y]+=value;
	double rand_value=RandomNumber(par.max_matrix);
	if (0.999<=rand_value/par.max_matrix){
		sigma[l][x][y]=1.;}
    // do not allow young lattice sites outside the cells
    if (dish->getCell(cpm->Sigma(x,y)).getTau() ==0 )
      sigma[l][x][y]= 0.;

    }

}

void PDE::MILayerCA(int l, double value, CellularPotts *cpm, Dish *dish){
  // First set pixels outside cells to 0, and new pixels within cells to 1
  {int new_sigma[sizex][sizey];
    for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++) {
    if (dish->getCell(cpm->Sigma(x,y)).getTau() ==0 && sigma[l][x][y]!=0){
    	  sigma[l][x][y]= 0;}
    else if (cpm->Sigma(x,y)>=1 && sigma[l][x][y]==0){
	sigma[l][x][y]=1;

}}

	//Do Eden growth
   for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++) {
	new_sigma[x][y]=sigma[l][x][y];
	int k;

	if (sigma[l][x][y]==1) {
	  // take a random neighbour
	  int xyp=(int)(8*RANDOM()+1);
	  int xp = nx[xyp]+x;
	  int yp = ny[xyp]+y;
	  int kp;

	//Check if both sites are part of the same cell
	if (dish->getCell(cpm->Sigma(x,y)).getTau()==dish->getCell(cpm->Sigma(xp,yp)).getTau()){
	  //  NB removing this border test yields interesting effects :-)
	  // You get a ragged border, which you may like!
		if ((kp=sigma[l][xp][yp])!=-1){
			if (kp>1){
			//Make site part of adhesion complex if next to adhesion complex
		// Probabililty depends on number of 'old' neighbours
		int old_neighbours=0;
     	for (int i1=-1;i1<=1;i1++)
     		for (int i2=-1;i2<=1;i2++){
			if (sigma[l][x+i1][y+i2]>1)
				old_neighbours+=1;}

				double random_double=rand()/double(RAND_MAX);
				if (random_double<par.eden_p){//old_neighbours){
	      			new_sigma[x][y]=2;}}
	    		else{
	     			new_sigma[x][y]=1;
}}
		else{
			//cout<< x <<", " << y << " border"<< endl;
			new_sigma[x][y]=1;}
	}
	//spontaneous formation of new adhesion
	double rand_spon = rand()/double(RAND_MAX);
	if (rand_spon < par.spontaneous_p){
		new_sigma[x][y]=2;}

}

	else if ((k=sigma[l][x][y])>1){//cpm->Sigma(x,y)>0){
		if (k<par.max_matrix-1){
		//cout << "1 " << new_sigma[x][y] << endl;
		new_sigma[x][y]=k+value;}
		//cout << "2 " << new_sigma[x][y] << endl;

		int young_neighbours=0;
     	for (int i1=-1;i1<=1;i1++)
     		for (int i2=-1;i2<=1;i2++){
			if (sigma[l][x+i1][y+i2]==1)
				young_neighbours+=1;}

	if (young_neighbours>0){	//cout<< x <<", " << y << " in cell old"<< endl;
		double rand_double=rand()/double(RAND_MAX);
		//cout << rand_double << endl;
		if (rand_double<par.decay_p){
			//cout<< x <<", " << y << " rejuvenate"<< endl;
			new_sigma[x][y]=1;}
		//else{
			//cout<< x <<", " << y << " stay the same"<< endl;
		//	new_sigma[x][y]= sigma[l][x][y]+1;}
	}}
      }
    // copy sigma to new_sigma, but do not touch the border!
	  {  for (int x=1;x<sizex-1;x++) {
      for (int y=1;y<sizey-1;y++) {
	//cout << x <<", " << y << "apply new sigmas" << endl;
	sigma[l][x][y]=new_sigma[x][y];
      }
    }
  }}



}
