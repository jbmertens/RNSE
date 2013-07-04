#ifndef RNSFLUID_H
#define RNSFLUID_H

/* evolution functions for different DOF's */
void jsource(PointData *paq);
void Jsource(PointData *paq);
static inline simType energy_evfn(PointData *paq);
static inline simType fluid_evfn(PointData *paq, int u);
static inline simType field_evfn(PointData *paq);
static inline simType ddtfield_evfn(PointData *paq);

void g2wevolve(simType *grid, simType *wedge, PointData *paq, int i, int j, int k);
void w2pevolve(simType *grid, simType *wedge, PointData *paq, int i, int j, int k);
void calculatequantities(PointData *paq, int i, int j, int k);


/* Deal with fluid coupling. */
static inline simType getXI();
static inline void setXI(simType xi);
static simType XI = 0;

/* get/set fluid coupling parameter. */
static inline void setXI(simType xi)
{
    XI = xi;
}
static inline simType getXI()
{
    return XI;
}


/* Deal with potential parameter. */
static inline simType getALPHA();
static inline void setALPHA(simType alpha);
static simType ALPHA = 0.90;

/* get/set potential parameter. */
static inline void setALPHA(simType alpha)
{
    ALPHA = alpha;
}
static inline simType getALPHA()
{
    return ALPHA;
}



#endif
