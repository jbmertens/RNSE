#ifndef IO_H
#define IO_H

/* io functionality prototypes */
typedef struct {
  char *data_dir;
  char *data_name;
  int fwrites;
  int datasize;
} IOData;

void writeinfo(IOData filedata);
void dumpstate(simType *fields, IOData filedata);
void readstate(simType *fields, IOData filedata);

#endif
