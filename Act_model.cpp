
/*

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
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
#include "info.h"
#include "parameter.h"
#include "sqr.h"
#include "leonie.cpp"



#include "qtgraph.h"
/*#else
#include "x11graph.h"
#endif*/



using namespace std;

void CellularPotts::ChangeThetas(Dish *dish);

INIT {

  try {

    // Define initial distribution of cells
    cout<< "Initialization"<< endl;
    // cout << "Pillar present " << CPM->AnyPillar() << endl;
    CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
// cout << "Pillar present " << CPM->AnyPillar() << endl;
// CPM->ReadZygotePicture();
    CPM->ConstructInitCells(*this, 1, par.target_area, par.target_perimeter);
// cout << "Pillar present " << CPM->AnyPillar() << endl;
// cout << "Pillar nearby (0,0)" << CPM->IsPillar(0,0)<< endl;
CPM->InitializeEdgeList();
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }

}

TIMESTEP {

  try {

    // cout << "new timestep" << endl;

    static int i=0;

    static Dish *dish=new Dish();
    static Info *info=new Info(*dish, *this);

    dish->CPM->AmoebaeMove(dish->PDEfield);

    dish->CPM->concentration_();

    dish->CPM->diffuse();

    if (par.max_Act){
      dish->PDEfield->AgeLayer(2,1.,dish->CPM, dish);
    }
    if (par.lambda_persistence){
      dish->CPM->ChangeThetas(dish);
    }

    char buff[400];
  	 snprintf(buff, sizeof(buff), "%s/cellcenter.txt",par.datadir);
  	 std::string buffAsStdStr = buff;
 	//cout<<buffAsStdStr<<"\n";
   	 std::ofstream out(buff, ios::app);
         info->WriteCOMsTorus(out);

    if (par.graphics && !(i%par.storage_stride)) {

      char fname[200];
      sprintf(fname,"%s/extend%05d.png",par.datadir,i);

      BeginScene();
      ClearImage();
      dish->CPM->PlotChem(this);


      dish->Plot(this);

      dish->PDEfield->PlotInCells (this, dish->CPM, 2);
      dish->CPM->SearchNandPlotClear(this);

      EndScene();
      info->Menu();


      Write(fname);

    }
    i++;
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }
}




int PDE::MapColour(double val) {

  return (((int)((val/((val)+1.))*100))%100)+155;
}


int PDE::MapColour3(double val, int l) {
	int step = (240)/par.max_Act;
  return (int)(256-val*step-1);
}


int main(int argc, char *argv[]) {


  try {

    QApplication a(argc, argv);

    // Read parameters
    par.Read(argv[1]);

    Seed(par.rseed);

    //QMainWindow mainwindow w;


    QtGraphics g(par.sizex*2,par.sizey*2);
		// cout<< "After  QtGraphics g(par.sizex*2,par.sizey*2);"<< endl;
    a.setMainWidget( &g );
    a.connect(&g, SIGNAL(SimulationDone(void)), SLOT(quit(void)) );

    if (par.graphics)
      g.show();

    a.exec();

  } catch(const char* error) {
    std::cerr << error << "\n";
    exit(1);
  }
  catch(...) {
    std::cerr << "An unknown exception was caught\n";
  }
  return 0;
}
