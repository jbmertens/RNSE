#include "defines.h"


/* 
 * Write simulation information to file.
 */
void writeinfo(IOData filedata)
{
  char *infofile;
  infofile = (char *) malloc(200 * sizeof(char));
  FILE *datafile;

  /* write simulation parameters to file */
  strcpy(infofile, filedata.data_dir);
  strcat(infofile, filedata.data_name);
  strcat(infofile, ".info");
  datafile = fopen(infofile, "w+");
  if( datafile == NULL )
  {
    printf("Error opening file: %s\n", infofile );
  }
  fprintf(datafile, "# Writing simulation.\n");
  fprintf(datafile, "# Points Sampled:\n%d\n", POINTS_TO_SAMPLE);
  fprintf(datafile, "# Steps recorded (requires (steps/recorded) to be an integer):\n%d\n",
    STEPS_TO_SAMPLE);
  fprintf(datafile, "# dx/dt (timestep):\n%f\n", dx/dt);
  fprintf(datafile, "# dx (voxel size):\n%f\n", dx);
  fprintf(datafile, "# Equation of state parameter w:\n%f\n", W_EOS);
  fprintf(datafile, "# Physical lattice dimensions:\n%f\t%f\t%f\n",
    SIZE, SIZE, SIZE);
  fprintf(datafile, "# Maximum number of steps run:\n%d\n", MAX_STEPS);
  fprintf(datafile, "# Initial bubble radius:\n%f\t%f\n",
    R0 /* physical size */, R0/SIZE /* size in pixels */);
  fprintf(datafile, "# Coupling Constant:\n%f\n", getXI());

  fclose(datafile);
  free(infofile);
}


/* 
 * Write current simulation state to file.  Should be able to handle arbitrary sized arrays.
 * The "datasize" parameter should be length of a cube.
 */
void dumpstate(simType *fields, IOData filedata)
{
  char *filename, *buffer;
  filename = (char *) malloc(100 * sizeof(char));
  buffer = (char *) malloc(5 * sizeof(char));

  /* file data for files */
  sprintf(buffer, "%d", filedata.fwrites);
  strcpy(filename, filedata.data_dir);
  strcat(filename, filedata.data_name);
  strcat(filename, ".");
  strcat(filename, buffer);
  strcat(filename, ".h5.gz");

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[RANK] = {filedata.datasize, filedata.datasize, filedata.datasize, DOF},
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

  return;
}


/* 
 * Dump a strip of data along an axis - output just scalar field values for now.
 */
void dumpstrip(simType *fields, IOData filedata)
{
  int i;

  // x-axis data
  char *filename, *buffer;
  gzFile *datafile;

  filename = malloc(200 * sizeof(*filename));
  strcpy(filename, filedata.data_dir);
  strcat(filename, filedata.data_name);
  strcat(filename,  ".strip.dat.gz");
  buffer = malloc(20 * sizeof(*buffer));

  datafile = (gzFile *)gzopen(filename, "ab");
  if(datafile == Z_NULL) {
    printf("Error opening file: %s\n", filename);
    return;
  }

  for(i=0; i<POINTS; i++)
  {
    // field values
    sprintf(buffer, "%g\t", fields[INDEX(i,POINTS/2,POINTS/2,1)]);
    gzwrite(datafile, buffer, strlen(buffer));
    // fluid velocity in direction of slice
    sprintf(buffer, "%g\t", fields[INDEX(POINTS/2,i,POINTS/2,2)]);
    gzwrite(datafile, buffer, strlen(buffer));
    // energy density
    sprintf(buffer, "%g\t", fields[INDEX(i,POINTS/2,POINTS/2,0)]);
    gzwrite(datafile, buffer, strlen(buffer));
  }  
  gzwrite(datafile, "\n", strlen("\n")); 
 
  gzclose(datafile);
  free(filename);
  free(buffer);

  return;
}


/*
 * Read in hdf5 field state.
 * Resets all other values in fields array.
 */
void readstate(simType *fields, IOData filedata)
{

  hid_t           file, dset, dcpl;    /* Handles */
  herr_t          status;
  htri_t          avail;
  H5Z_filter_t    filter_type;

  size_t          nelmts;
  unsigned int    flags,
                  filter_info;
  int             i, j, k, l;

  simType         min;

  /*
   * Open file and dataset using the default properties.
   */
  printf("Attempting to open file: %s\n", filedata.read_data_name);

  file = H5Fopen(filedata.read_data_name, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset = H5Dopen2(file, "Dataset1", H5P_DEFAULT);
  dcpl = H5Dget_create_plist(dset);

  /*
   * Read data in and find the maximum value in the dataset, to
   * verify that it was read correctly.
   */
  status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fields);
  printf ("Status after reading is: %d\n", status);

  /*
   * Close and release resources.
   */
  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Fclose (file);

  // Expand data to fill array appropriately (work backwords so data isn't overwritten)
  for(i=POINTS-1; i>=0; i--) {
    for(j=POINTS-1; j>=0; j--)
      for(k=POINTS-1; k>=0; k--)
      {
        fields[INDEX(i,j,k,4)] = fields[POINTS*POINTS*(i) + POINTS*(j) + (k)];
      }
    }
  LOOP3(i,j,k)
  {
    fields[INDEX(i,j,k,0)] = LOG_E;
    fields[INDEX(i,j,k,1)] = 0.0;
    fields[INDEX(i,j,k,2)] = 0.0;
    fields[INDEX(i,j,k,3)] = 0.0;
    fields[INDEX(i,j,k,5)] = 0.0;
  }

  return;
}

/* 
 * Explicitly write what physical time things occurred at - allow variable timestep.
 */
void write_timestep(double deltat, int stepnum, IOData filedata)
{
  char *infofile;
  infofile = (char *) malloc(200 * sizeof(char));
  FILE *datafile;

  /* write simulation parameters to file */
  strcpy(infofile, filedata.data_dir);
  strcat(infofile, filedata.data_name);
  strcat(infofile, ".timeinfo");
  datafile = fopen(infofile, "a+");
  if( datafile == NULL )
  {
    printf("Error opening file: %s\n", infofile );
  }

  /* write data here! */
  fprintf(datafile, "%d\t%10.10f\n", stepnum, deltat);

  fclose(datafile);
  free(infofile);
}
