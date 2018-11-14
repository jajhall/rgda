#ifndef HFACTOR_H_
#define HFACTOR_H_

#include "HMatrix.h"
#include "HVector.h"

#include <vector>
using namespace std;

class HFactor {
public:
    void setup(const HMatrix *matrix);
    void build(int *Bindex);
    void ftran(HVector& vector, double hist_dsty) const;
    void btran(HVector& vector, double hist_dsty) const;
    void update(HVector& column, HVector& row_ep, int iRow, bool& re_inv);

private:
    void build_prepare(const int *base);
    void build_singleton();
    void build_singleton_column();
    void build_singleton_row();
    void build_kernel_prepare();
    bool build_kernel_search();
    void build_kernel_eliminate_prepare();
    void build_kernel_eliminate_columns();
    void build_kernel_eliminate_rows();
    void build_expand_column(const int iCol, const int fillin);
    void build_expand_row(const int iRow);
    void build_finish(int *base);

    void ftran_L(HVector& vector) const;
    void ftran_R(HVector& vector) const;
    void ftran_U(HVector& vector) const;
    void btran_U(HVector& vector) const;
    void btran_R(HVector& vector) const;
    void btran_L(HVector& vector) const;

    void ftran_LH(HVector& vector) const;
    void ftran_UH(HVector& vector) const;
    void btran_UH(HVector& vector) const;
    void btran_LH(HVector& vector) const;

    // Search parameters
    double search_tol;
    int search_lim;

    // Problem size and (pointer to) coefficient matrix
    int numRow, numCol;
    const HMatrix *matrix;

    // Working buffer
    int nwork;
    vector<int> iwork_;

    // Basis matrix, Factor L and U, update matrix
    int BlimitX, BnX, LnX, UnX, PFnX;
    vector<int> Bstart_, Lstart_, Ustart_, LRstart_, URstart_, PFstart_;
    vector<int> Bindex_, Lindex_, Uindex_, LRindex_, URindex_, PFindex_;
    vector<double> Bvalue_, Lvalue_, Uvalue_, LRvalue_, URvalue_, PFvalue_;

    // All pivots records
    int pivotK, pivotCol, pivotRow, npivot, PFnPivot;
    vector<int> permute_, HpivotIndex_;
    vector<double> HpivotValue_, PFpivotValue_;

    // Markowitz kernel matrices
    int MCnX, MRnX, MRnFillin;
    HVector Mcolumn;
    vector<int> MCstart, MCcountA, MCplaceN, MCspace;
    vector<int> MRstart, MRcount, MRspace, MRcountb4;
    vector<int> MCindex, MRindex, MRfillin;
    vector<double> MCminpivot, MCvalue;

    // Count-link-list
    vector<int> clinknext, clinklast, clink1;
    vector<int> rlinknext, rlinklast, rlink1;

    // Operations on the count-link-list
    void clinkAdd(const int index, const int count) {
        const int mover = clink1[count];
        clinklast[index] = -2 - count;
        clinknext[index] = mover;
        clink1[count] = index;
        if (mover >= 0)
            clinklast[mover] = index;
    }
    void clinkDel(const int index) {
        const int xlast = clinklast[index];
        const int xnext = clinknext[index];
        if (xlast >= 0)
            clinknext[xlast] = xnext;
        else
            clink1[-xlast - 2] = xnext;
        if (xnext >= 0)
            clinklast[xnext] = xlast;
    }

    void rlinkAdd(const int index, const int count) {
        const int mover = rlink1[count];
        rlinklast[index] = -2 - count;
        rlinknext[index] = mover;
        rlink1[count] = index;
        if (mover >= 0)
            rlinklast[mover] = index;
    }
    void rlinkDel(const int index) {
        const int xlast = rlinklast[index];
        const int xnext = rlinknext[index];
        if (xlast >= 0)
            rlinknext[xlast] = xnext;
        else
            rlink1[-xlast - 2] = xnext;
        if (xnext >= 0)
            rlinklast[xnext] = xlast;
    }
};

#endif /* HFACTOR_H_ */
