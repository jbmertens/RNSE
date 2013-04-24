#ifndef IO_H
#define IO_H

/* io functionality prototypes */
void writeinfo(IOData filedata);
void dumpstate(simType *fields, IOData filedata);
void dumpstrip(simType *fields, IOData filedata);
void readstate(simType *fields, IOData filedata);
void write_timestep(double deltat, int stepnum, IOData filedata);

#endif
