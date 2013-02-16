
#include "defines.h"

/* 
 * Write simulation information to file.
 */
void writeinfo(char *dir, char *root)
{
  char *infofile;
  infofile = (char *) malloc(200 * sizeof(char));
  FILE *datafile;

  /* write simulation parameters to file */
  strcpy(infofile, dir);
  strcat(infofile, root);
  strcat(infofile, ".info");
  datafile = fopen(infofile, "w+");
  if( datafile == NULL )
  {
    printf("Error opening file: %s\n", infofile );
  }
  fprintf(datafile, "# Writing simulation.\n");
  fprintf(datafile, "#Points Sampled:\n%d\n", POINTS_TO_SAMPLE);
  fprintf(datafile, "# Steps recorded:\n%d\n", STEPS_TO_RECORD);
  fprintf(datafile, "# dx/dt (timestep):\n%f", dx/dt);
  fprintf(datafile, "# Equation of state parameter w:\n%f", W_EOS);
  fprintf(datafile, "# Physical lattice dimensions:\n%f %f %f\n",
    SIZE, SIZE, SIZE);
  fprintf(datafile, "# Total number of steps run:\n%d\n", STEPS);

  fclose(datafile);
  free(infofile);
}


/* 
 * Write current simulation state to file.  Should be able to handle arbitrary sized arrays.
 * The "datasize" parameter should be length of a cube.
 */
void dumpstate(simType *fields, int fwrites, int datasize, char *dir, char *name)
{
  char *filename, *buffer;
  filename = (char *) malloc(100 * sizeof(char));
  buffer = (char *) malloc(5 * sizeof(char));

  /* file data for files */
  sprintf(buffer, "%d", fwrites);
  strcpy(filename, dir);
  strcat(filename, name);
  strcat(filename, ".");
  strcat(filename, buffer);
  strcat(filename, ".h5.gz");

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[RANK] = {datasize, datasize, datasize, DOF},
          maxdims[RANK] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED},
          chunk[RANK] = {6, 6, 6, 6};

  file = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (RANK, dims, maxdims);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_deflate (dcpl, 9);
  status = H5Pset_chunk (dcpl, RANK, chunk);
  dset = H5Dcreate2 (file, "DS1", H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fields);

  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

  free(filename);
  free(buffer);
}

// need to write this function.
void readstate(simType *fields, int fwrites)
{

}
