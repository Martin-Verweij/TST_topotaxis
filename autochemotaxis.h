#include "parameter.h"

float concentration[500][500] = { 0 };

void concentration_()
{
  int i, j,q;

  for ( i = 0; i < sizex-1; i++ )
    for ( j = 0; j < sizey-1; j++ ) {

      if (sigma[i][j]<=0) {      }

      else {

      //  C.concentration[i][j] += par.secr_rate_;

      }

      if ( sigma[i][j] != sigma[i+1][j] ) {

        concentration[i][j] += par.secr_rate_;

      }

      if ( sigma[i][j] != sigma[i][j+1] ) {

        concentration[i][j] += par.secr_rate_;

      }

      if (sigma[i][j]!=sigma[i+1][j+1] || sigma[i+1][j]!=sigma[i][j+1] ) {

        concentration[i][j] += par.secr_rate_;

      }
      // else
      //   if (g && sigma[i][j]>0)
      //     g->Point( colour, 2*i+1, 2*j+1 );
    }

}
