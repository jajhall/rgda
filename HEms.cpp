#include "HModel.h"
#include "HConst.h"
#include <cstdio>
#include <cmath>

HModel hmodel;

const double HEMS_INF = 1e31;

#ifdef __GNUC__
#define HEMS_API void
#else
#define HEMS_API extern "C" __declspec(dllexport) void __stdcall
#endif

HEMS_API hems_dsca(int *rtcod, double *dspace, int ndwords, int nmodel) {
    if (hmodel.solving_with_log)
        puts("hems_dsca");

    *rtcod = 0;
}

HEMS_API hems_init(int *rtcod, double *dspace) {
    if (hmodel.solving_with_log)
        puts("hems_init");

    *rtcod = 0;
}

HEMS_API hems_mset(int *rtcod,
                   double *dspace,
                   int strtnum,
                   int maxalw,
                   int maxprt,
                   int trace,
                   int userexit,
                   int endnum,
                   int nonum) {
    if (hmodel.solving_with_log)
        puts("hems_mset");

    *rtcod = 0;
}

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
                   double *Avalue) {
    if (hmodel.solving_with_log) {
        puts("hems_xmps");
        printf("    Allocated size = %6d columns %6d rows and %8d elements\n",
               *numCol, *numRow, *AcountX);
    }

    HInput hinput;
    hinput.loadMps("test.mps");
    int numCol_ = hinput.numCol;
    int numRow_ = hinput.numRow;
    int AcountX_ = hinput.AcountX;

    if (hmodel.solving_with_log) {
        printf("    MPS model size = %6d columns %6d rows and %8d elements\n",
               numCol_, numRow_, AcountX_);
    }

    if (*numRow < numRow_ || *numCol < numCol_ || *AcountX < AcountX_) {
        // There is not enough spaces
        *rtcod = 1;
        return;
    }

    *numRow = numRow_;
    *numCol = numCol_;
    *AcountX = AcountX_;

    // Copy matrix (while change to FORTRAN convention)
    for (int i = 0; i <= numCol_; i++) {
        Astart[i] = hinput.Astart[i] + 1;
    }
    for (int i = 0; i < AcountX_; i++) {
        Aindex[i] = hinput.Aindex[i] + 1;
        Avalue[i] = hinput.Avalue[i];
    }

    // Copy cost and bounds (while change to EMSOL infinity)
    for (int i = 0; i < numCol_; i++) {
        cost[i] = hinput.cost[i];
        colLower[i] = max(hinput.colLower[i], -HEMS_INF);
        colUpper[i] = min(hinput.colUpper[i], +HEMS_INF);
    }
    for (int i = 0; i < numRow_; i++) {
        rowLower[i] = max(hinput.rowLower[i], -HEMS_INF);
        rowUpper[i] = min(hinput.rowUpper[i], +HEMS_INF);
    }

    // Exit successfully
    *rtcod = 0;

}

#include <iostream>
using namespace std;

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
                   double *Avalue) {
    if (hmodel.solving_with_log)
        puts("hems_lmdl");

    HInput& hinput = hmodel.hinput;
    hinput.numCol = numCol;
    hinput.numRow = numRow;
    hinput.AcountX = AcountX;

    hinput.Astart.assign(Astart, Astart + numCol + 1);
    hinput.Aindex.assign(Aindex, Aindex + AcountX);
    hinput.Avalue.assign(Avalue, Avalue + AcountX);
    for (int i = 0; i <= numCol; i++)
        hinput.Astart[i]--;
    for (int k = 0; k < AcountX; k++)
        hinput.Aindex[k]--;

    // Remove explicit elements
    int kput = 0;
    for (int i = 0; i < numCol; i++) {
        int start = hinput.Astart[i];
        int end = hinput.Astart[i + 1];
        // Update last pointer
        if (start > kput)
            hinput.Astart[i] = kput;
        for (int k = start; k < end; k++) {
            if (fabs(hinput.Avalue[k]) > H_TT) {
                hinput.Aindex[kput] = hinput.Aindex[k];
                hinput.Avalue[kput] = hinput.Avalue[k];
                kput++;
            }
        }
    }
    // Update the matrix size
    if (kput != AcountX) {
        hinput.AcountX = kput;
        hinput.Astart[numCol] = kput;
        hinput.Aindex.resize(kput);
        hinput.Avalue.resize(kput);
    }

    // Copy bounds
    hinput.colLower.assign(colLower, colLower + numCol);
    hinput.colUpper.assign(colUpper, colUpper + numCol);
    hinput.rowLower.assign(rowLower, rowLower + numRow);
    hinput.rowUpper.assign(rowUpper, rowUpper + numRow);
    hinput.cost.assign(cost, cost + numCol);

    // Change infinity to HSOL infinity
    for (int i = 0; i < numCol; i++) {
        if (hinput.colLower[i] < -HEMS_INF * 0.99)
            hinput.colLower[i] = -H_INF;
        if (hinput.colUpper[i] > +HEMS_INF * 0.99)
            hinput.colUpper[i] = +H_INF;
    }
    for (int i = 0; i < numRow; i++) {
        if (hinput.rowLower[i] < -HEMS_INF * 0.99)
            hinput.rowLower[i] = -H_INF;
        if (hinput.rowUpper[i] > +HEMS_INF * 0.99)
            hinput.rowUpper[i] = +H_INF;
    }

    // Setup other part of the model
    hmodel.setup();
    *rtcod = 0;
}

HEMS_API hems_scal(int *rtcod, double *dspace) {
    if (hmodel.solving_with_log)
        puts("hems_scal");
    hmodel.scale();
    *rtcod = 0;
}

HEMS_API hems_sslv(int *rtcod, double *dspace, int algo, int init) {
    if (hmodel.solving_with_log)
        puts("hems_sslv");
    hmodel.solve();
    *rtcod = 0;
}

HEMS_API hems_iget(int *rtcod, double *dspace, int *iarray, int number) {
    if (hmodel.solving_with_log)
        puts("hems_iget");
    iarray[46] = hmodel.problem_status;
    *rtcod = 0;
}

HEMS_API hems_rget(int *rtcod, double *dspace, double *darray, int number) {
    if (hmodel.solving_with_log)
        puts("hems_rget");
    darray[17] = hmodel.objective;
    *rtcod = 0;
}

HEMS_API hems_gepr(int *rtcod,
                   double *dspace,
                   int rowListSize,
                   int colListSize,
                   int *rowList,
                   int *colList,
                   double *rowValue,
                   double *colValue) {

    if (hmodel.solving_with_log)
        puts("hems_gepr");

    int numRow = hmodel.houtput.numRow;
    int numCol = hmodel.houtput.numCol;
    double *value = &hmodel.houtput.value[0];
    if (rowListSize > numRow) {
        for (int i = 0; i < numRow; i++)
            rowValue[i] = -value[i + numCol];
    } else {
        for (int i = 0; i < rowListSize; i++)
            rowValue[i] = -value[rowList[i] + numCol - 1];
    }
    if (colListSize > numCol) {
        for (int i = 0; i < numCol; i++)
            colValue[i] = value[i];
    } else {
        for (int i = 0; i < colListSize; i++)
            colValue[i] = value[colList[i] - 1];
    }

    *rtcod = 0;
}

HEMS_API hems_gedu(int *rtcod,
                   double *dspace,
                   int rowListSize,
                   int colListSize,
                   int *rowList,
                   int *colList,
                   double *rowDual,
                   double *colDual) {
    if (hmodel.solving_with_log)
        puts("hems_gedu");

    int numRow = hmodel.houtput.numRow;
    int numCol = hmodel.houtput.numCol;
    double *dual = &hmodel.houtput.dual[0];

    if (rowListSize > numRow) {
        for (int i = 0; i < numRow; i++)
            rowDual[i] = dual[i + numCol];
    } else {
        for (int i = 0; i < rowListSize; i++)
            rowDual[i] = dual[rowList[i] + numCol - 1];
    }

    if (colListSize > numCol) {
        for (int i = 0; i < numCol; i++)
            colDual[i] = dual[i];
    } else {
        for (int i = 0; i < colListSize; i++)
            colDual[i] = dual[colList[i] - 1];
    }
    *rtcod = 0;
}

HEMS_API hems_rgda(int *rtcod, double *dspace) {
    if (hmodel.solving_with_log)
        puts("hems_rgda");

    hmodel.sense();
    *rtcod = 0;
}

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
                   int *rb_dn_l) {

    if (hmodel.solving_with_log)
        puts("hems_gerg");

    HOutput *output = &hmodel.houtput;
    int numRow = output->numRow;
    int numCol = output->numCol;
    if (rowListSize > numRow) {
        for (int i = 0; i < numRow; i++) {
            int iTot = i + numCol;
            // Copy all row bounds
            rb_up_b[i] = -output->b_dn_b[iTot];
            rb_up_f[i] = output->b_dn_f[iTot];
            rb_up_e[i] = output->b_dn_e[iTot] + 1;
            rb_up_l[i] = output->b_dn_l[iTot] + 1;
            rb_dn_b[i] = -output->b_up_b[iTot];
            rb_dn_f[i] = output->b_up_f[iTot];
            rb_dn_e[i] = output->b_up_e[iTot] + 1;
            rb_dn_l[i] = output->b_up_l[iTot] + 1;
        }
    } else {
        for (int i = 0; i < rowListSize; i++) {
            int iRow = rowList[i] - 1;
            int iTot = iRow + numCol;
            // Copy all row bounds
            rb_up_b[i] = -output->b_dn_b[iTot];
            rb_up_f[i] = output->b_dn_f[iTot];
            rb_up_e[i] = output->b_dn_e[iTot] + 1;
            rb_up_l[i] = output->b_dn_l[iTot] + 1;
            rb_dn_b[i] = -output->b_up_b[iTot];
            rb_dn_f[i] = output->b_up_f[iTot];
            rb_dn_e[i] = output->b_up_e[iTot] + 1;
            rb_dn_l[i] = output->b_up_l[iTot] + 1;
        }
    }

    // Change row bounds value to HEMS style

    if (colListSize > numCol) {
        for (int i = 0; i < numCol; i++) {
            // Copy all column SBND
            cb_up_b[i] = output->b_up_b[i];
            cb_up_f[i] = output->b_up_f[i];
            cb_up_e[i] = output->b_up_e[i] + 1;
            cb_up_l[i] = output->b_up_l[i] + 1;
            cb_dn_b[i] = output->b_dn_b[i];
            cb_dn_f[i] = output->b_dn_f[i];
            cb_dn_e[i] = output->b_dn_e[i] + 1;
            cb_dn_l[i] = output->b_dn_l[i] + 1;

            // Copy all column SOBJ
            co_up_c[i] = output->c_up_c[i];
            co_up_f[i] = output->c_up_f[i];
            co_up_e[i] = output->c_up_e[i] + 1;
            co_up_l[i] = output->c_up_l[i] + 1;
            co_dn_c[i] = output->c_dn_c[i];
            co_dn_f[i] = output->c_dn_f[i];
            co_dn_e[i] = output->c_dn_e[i] + 1;
            co_dn_l[i] = output->c_dn_l[i] + 1;
        }
    } else {
        for (int i = 0; i < colListSize; i++) {
            int iCol = colList[i] - 1;

            // Copy all column SBND
            cb_up_b[i] = output->b_up_b[iCol];
            cb_up_f[i] = output->b_up_f[iCol];
            cb_up_e[i] = output->b_up_e[iCol] + 1;
            cb_up_l[i] = output->b_up_l[iCol] + 1;
            cb_dn_b[i] = output->b_dn_b[iCol];
            cb_dn_f[i] = output->b_dn_f[iCol];
            cb_dn_e[i] = output->b_dn_e[iCol] + 1;
            cb_dn_l[i] = output->b_dn_l[iCol] + 1;

            // Copy all column SOBJ
            co_up_c[i] = output->c_up_c[iCol];
            co_up_f[i] = output->c_up_f[iCol];
            co_up_e[i] = output->c_up_e[iCol] + 1;
            co_up_l[i] = output->c_up_l[iCol] + 1;
            co_dn_c[i] = output->c_dn_c[iCol];
            co_dn_f[i] = output->c_dn_f[iCol];
            co_dn_e[i] = output->c_dn_e[iCol] + 1;
            co_dn_l[i] = output->c_dn_l[iCol] + 1;
        }
    }

    // Change to EMS convention
    int c_size = colListSize <= numCol ? colListSize : numCol;
    int r_size = rowListSize <= numRow ? rowListSize : numRow;
    double buffer_size[12] = { c_size, c_size, c_size, c_size, c_size, c_size,
            c_size, c_size, r_size, r_size, r_size, r_size };
    double *dbuffer[12] = { co_up_c, co_dn_c, co_up_f, co_dn_f, cb_up_b,
            cb_dn_b, cb_up_f, cb_dn_f, rb_up_b, rb_dn_b, rb_up_f, rb_dn_f };
    int *ibuffer[12] = { co_up_e, co_dn_e, co_up_l, co_dn_l, cb_up_e, cb_dn_e,
            cb_up_l, cb_dn_l, rb_up_e, rb_dn_e, rb_up_l, rb_dn_l };

    // Change all value to EMS infinity style
    for (int i = 0; i < 12; i++) {
        double *dwork = dbuffer[i];
        int size = buffer_size[i];
        for (int k = 0; k < size; k++) {
            dwork[k] = max(dwork[k], -HEMS_INF);
            dwork[k] = min(dwork[k], +HEMS_INF);
        }
    }

    // Change all index to output requirement (iRow in [-numRow, -1])
    for (int i = 0; i < 12; i++) {
        int *iwork = ibuffer[i];
        int size = buffer_size[i];
        for (int k = 0; k < size; k++) {
            if (iwork[k] > numCol)
                iwork[k] = numCol - iwork[k];
        }
    }

    *rtcod = 0;
}
