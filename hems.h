#ifndef HEMS_H_
#define HEMS_H_

#ifdef __GNUC__
#define HEMS_API void
#else
#define HEMS_API extern "C" void __stdcall
#endif


HEMS_API hems_dsca(int *rtcod, double *dspace, int ndwords, int nmodel);

HEMS_API hems_init(int *rtcod, double *dspace);

HEMS_API hems_mset(int *rtcod,
                   double *dspace,
                   int strtnum,
                   int maxalw,
                   int maxprt,
                   int trace,
                   int userexit,
                   int endnum,
                   int nonum);

HEMS_API hems_xmps(int *rtcod,
                   double *dspace,
                   int type,
                   int *numRow,
                   int *numCol,
                   int *AcountX,
                   double *cost,
                   double *rowLower,
                   double *rowUpper,
                   double *colLower,
                   double *colUpper,
                   int *Aindex,
                   int *Astart,
                   double *Avalue);

HEMS_API hems_lmdl(int *rtcod,
                   double *dspace,
                   int type,
                   int numRow,
                   int numCol,
                   int AcountX,
                   double *cost,
                   double *rowLower,
                   double *rowUpper,
                   double *colLower,
                   double *colUpper,
                   int *Aindex,
                   int *Astart,
                   double *Avalue);

HEMS_API hems_scal(int *rtcod, double *dspace);

HEMS_API hems_sslv(int *rtcod, double *dspace, int algo, int init);

HEMS_API hems_iget(int *rtcod, double *dspace, int *iarray, int number);

HEMS_API hems_rget(int *rtcod, double *dspace, double *darray, int number) ;

HEMS_API hems_gepr(int *rtcod,
                   double *dspace,
                   int rowListSize,
                   int colListSize,
                   int *rowList,
                   int *colList,
                   double *rowValue,
                   double *colValue);

HEMS_API hems_gedu(int *rtcod,
                   double *dspace,
                   int rowListSize,
                   int colListSize,
                   int *rowList,
                   int *colList,
                   double *rowDual,
                   double *colDual);

HEMS_API hems_rgda(int *rtcod, double *dspace);

HEMS_API hems_gerg(int *rtcod,
                   double *dspace,
                   int rowListSize,
                   int colListSize,
                   int *rowList,
                   int *colList,
                   double *co_up_c,
                   double *co_dn_c,
                   double *co_up_f,
                   double *co_dn_f,
                   int *co_up_e,
                   int *co_dn_e,
                   int *co_up_l,
                   int *co_dn_l,
                   double *cb_up_b,
                   double *cb_dn_b,
                   double *cb_up_f,
                   double *cb_dn_f,
                   int *cb_up_e,
                   int *cb_dn_e,
                   int *cb_up_l,
                   int *cb_dn_l,
                   double *rb_up_b,
                   double *rb_dn_b,
                   double *rb_up_f,
                   double *rb_dn_f,
                   int *rb_up_e,
                   int *rb_dn_e,
                   int *rb_up_l,
                   int *rb_dn_l);

#endif

