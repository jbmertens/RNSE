#ifndef RNSFLUID_H
#define RNSFLUID_H

/* evolution functions for different DOF's */
void jsource(PointData *paq);
void Jsource(PointData *paq);
static inline simType energy_evfn(PointData *paq);
static inline simType fluid_evfn(PointData *paq, int u);
static inline simType field_evfn(PointData *paq);
static inline simType ddtfield_evfn(PointData *paq);
void calculatequantities(simType *fields, PointData *paq, int i, int j, int k);
void evolve(simType *initial, simType *final, simType coeff, PointData *paq, int i, int j, int k);

/* Deal with fluid coupling. */
static inline simType getXI();
static inline void setXI(simType xi);
static simType XI;

/* get/set fluid coupling parameter. */
static inline void setXI(simType xi)
{
    XI = xi;
}
static inline simType getXI()
{
    return XI;
}

#endif
