#include "defines.h"


void hartleydump(simType *fields, simType *storage, IOData filedata)
{
  /* set values in storage array - storage should still contain STORAGE amount of storage. */
  int i, j, k, d;
  for(i=0; i<POINTS; i++)
    for(j=0; j<POINTS; j++)
      for(k=0; k<POINTS; k++)
        storage[POINTS*POINTS*i + POINTS*j + k] = fields[INDEX(i, j, k, 0)];

  /* create plan */
  fftw_plan p;
  p = fftw_plan_r2r_3d(POINTS, POINTS, POINTS,
                            storage, storage,
                            FFTW_DHT, FFTW_DHT, FFTW_DHT,
                            FFTW_ESTIMATE|FFTW_DESTROY_INPUT);

  /* Execute! */
  fftw_execute(p);

  /* Storage array should now contain DHT values. To save, we can rearrange values, creating a copy
     of fields, but with the field derivative replaced by DHT values.  By working from the end back
     in, we don't blow away any DHT values.  */
  for(i=POINTS-1; i>=0; i--)
    for(j=POINTS-1; j>=0; j--)
      for(k=POINTS-1; k>=0; k--)
        for(d=DOF-1; d>=0; d--)
        {
          if(5 == d) {
            storage[INDEX(i, j, k, 5)] = storage[POINTS*POINTS*i + POINTS*j + k];
          }
          else
          {
            storage[INDEX(i, j, k, d)] = fields[INDEX(i, j, k, d)];
          }
        }

  IOData dumpfile;
  dumpfile.data_dir = filedata.data_dir;
  dumpfile.data_name = "DUMP";
  dumpfile.fwrites = filedata.fwrites;
  dumpfile.datasize = POINTS;
  dumpstate(storage, dumpfile);

  return;
}
