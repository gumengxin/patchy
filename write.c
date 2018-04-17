/***************************************
 * write movie snapshot in .xyz format
 *
 * number of particles
 * blank
 * name x y z
 * C 0.0000 0.0000 0.0000
 * S 0.0000 1.0000 0.2345
 *
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include "system.h"

void Writemovie(FILE *FilePtr)
{
  int i,j;
  int pbccheckx,pbcchecky,pbccheckz;
  double rx,ry,rz;

  fprintf(FilePtr,"%d\n",NA+NB+NA*ZA+NB*ZB);
  fprintf(FilePtr,"%lf %lf %lf %lf %lf\n",LX,LY,LZ,Tinstant,Pinstant);

  for(i=0;i<NumberOfParticles;i++)
  {
    if(identity[i] == 1)	
     fprintf(FilePtr,"%s\t","N");
    else	
     fprintf(FilePtr,"%s\t","S");

    rx = PBC(RX[i],LX);
    ry = PBC(RY[i],LY);
    rz = PBC(RZ[i],LZ);

    if(rx != RX[i]) pbccheckx = 1;
    else pbccheckx = 0;
    if(ry != RY[i]) pbcchecky = 1;
    else pbcchecky = 0;
    if(rz != RZ[i]) pbccheckz = 1;
    else pbccheckz = 0;

    fprintf(FilePtr,"%lf\t",PBC(RX[i],LX));
    fprintf(FilePtr,"%lf\t",PBC(RY[i],LY));
    fprintf(FilePtr,"%lf\n",PBC(RZ[i],LZ));

   for(j=0;j<zcoordinate[i];j++)
   {
    if(identity[i] == 1)	
    fprintf(FilePtr,"%s\t","H");
    else
    fprintf(FilePtr,"%s\t","O");

    if(pbccheckx == 1)
    fprintf(FilePtr,"%lf\t",RSX[i][j]+rx-RX[i]);
    else
    fprintf(FilePtr,"%lf\t",RSX[i][j]);
    if(pbcchecky == 1)
    fprintf(FilePtr,"%lf\t",RSY[i][j]+ry-RY[i]);
    else
    fprintf(FilePtr,"%lf\t",RSY[i][j]);
    if(pbccheckz == 1)
    fprintf(FilePtr,"%lf\n",RSZ[i][j]+rz-RZ[i]);
    else
    fprintf(FilePtr,"%lf\n",RSZ[i][j]);
   }
    
  }

 // for(i=0;i<NumberOfParticles;i++)

}
