#include "HMatrix.h"
#include "HConst.h"
#include <cmath>

void HMatrix::setup(int numCol,
                    int numRow,
                    const int *Astart,
                    const int *Aindex,
                    const double *Avalue) {
    // Copy A
    numCol_ = numCol;
    numRow_ = numRow;
    Astart_.assign(Astart, Astart + numCol + 1);

    int AcountX = Astart[numCol];
    Aindex_.assign(Aindex, Aindex + AcountX);
    Avalue_.assign(Avalue, Avalue + AcountX);

    // Build row copy - pointers
    ARstart_.resize(numRow + 1);
    ARend_n_.assign(numRow, 0);
    for (int k = 0; k < AcountX; k++)
        ARend_n_[Aindex_[k]]++;
    ARstart_[0] = 0;
    for (int i = 1; i <= numRow_; i++)
        ARstart_[i] = ARstart_[i - 1] + ARend_n_[i - 1];
    for (int i = 0; i < numRow_; i++)
        ARend_n_[i] = ARstart_[i];

    // Build row copy - elements
    ARindex_.resize(AcountX);
    ARvalue_.resize(AcountX);
    for (int iCol = 0; iCol < numCol_; iCol++) {
        for (int k = Astart_[iCol]; k < Astart_[iCol + 1]; k++) {
            int iRow = Aindex_[k];
            int iPut = ARend_n_[iRow]++;
            ARindex_[iPut] = iCol;
            ARvalue_[iPut] = Avalue_[k];
        }
    }
}

void HMatrix::scale(double *col_scale, double *row_scale) {
    for (int iCol = 0; iCol < numCol_; iCol++)
        for (int k = Astart_[iCol]; k < Astart_[iCol + 1]; k++)
            Avalue_[k] *= (col_scale[iCol] * row_scale[Aindex_[k]]);

    for (int iRow = 0; iRow < numRow_; iRow++)
        for (int k = ARstart_[iRow]; k < ARstart_[iRow + 1]; k++)
            ARvalue_[k] *= (col_scale[ARindex_[k]] * row_scale[iRow]);

}

void HMatrix::update(int columnIn, int columnOut) {
    // 1. Move out some from basic part
    if (columnOut < numCol_) {
        for (int k = Astart_[columnOut]; k < Astart_[columnOut + 1]; k++) {
            int iRow = Aindex_[k];
            int iFind = ARend_n_[iRow];
            int iSwap = ARend_n_[iRow]++;
            while (ARindex_[iFind] != columnOut)
                iFind++;
            swap(ARindex_[iFind], ARindex_[iSwap]);
            swap(ARvalue_[iFind], ARvalue_[iSwap]);
        }
    }

    // 2. Move in some to the basic part
    if (columnIn < numCol_) {
        for (int k = Astart_[columnIn]; k < Astart_[columnIn + 1]; k++) {
            int iRow = Aindex_[k];
            int iFind = ARstart_[iRow];
            int iSwap = --ARend_n_[iRow];
            while (ARindex_[iFind] != columnIn)
                iFind++;
            swap(ARindex_[iFind], ARindex_[iSwap]);
            swap(ARvalue_[iFind], ARvalue_[iSwap]);
        }
    }
}

double HMatrix::compute_dot(HVector& vector, int iCol) {
    double result = 0;
    if (iCol < numCol_) {
        for (int k = Astart_[iCol]; k < Astart_[iCol + 1]; k++)
            result += vector.array[Aindex_[k]] * Avalue_[k];
    } else {
        result = vector.array[iCol - numCol_];
    }
    return result;
}

void HMatrix::collect_aj(HVector& vector, int iCol, double multi) {
    if (iCol < numCol_) {
        for (int k = Astart_[iCol]; k < Astart_[iCol + 1]; k++) {
            int index = Aindex_[k];
            double value0 = vector.array[index];
            double value1 = value0 + multi * Avalue_[k];
            if (value0 == 0)
                vector.index[vector.count++] = index;
            vector.array[index] = (fabs(value1) < H_TT) ? H_TZ : value1;
        }
    } else {
        int index = iCol - numCol_;
        double value0 = vector.array[index];
        double value1 = value0 + multi;
        if (value0 == 0)
            vector.index[vector.count++] = index;
        vector.array[index] = (fabs(value1) < H_TT) ? H_TZ : value1;
    }
}

void HMatrix::price_by_col(HVector& row_ap, HVector& row_ep) {
    // Alias
    int ap_count = 0;
    int *ap_index = &row_ap.index[0];
    double *ap_array = &row_ap.array[0];
    const double *ep_array = &row_ep.array[0];

    const int *Astart = &Astart_[0];
    const int *Aindex = &Aindex_[0];
    const double *Avalue = &Avalue_[0];

    // Computation
    for (int iCol = 0; iCol < numCol_; iCol++) {
        double value = 0;
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
            value += ep_array[Aindex[k]] * Avalue[k];
        }
        if (fabs(value) > H_TT) {
            ap_array[iCol] = value;
            ap_index[ap_count++] = iCol;
        }
    }
    row_ap.count = ap_count;

}

void HMatrix::price_by_row(HVector& row_ap, HVector& row_ep) {
    // Alias
    int ap_count = 0;
    int *ap_index = &row_ap.index[0];
    double *ap_array = &row_ap.array[0];
    const int ep_count = row_ep.count;
    const int *ep_index = &row_ep.index[0];
    const double *ep_array = &row_ep.array[0];

    const int *ARstart = &ARstart_[0];
    const int *ARend_n = &ARend_n_[0];
    const int *ARindex = &ARindex_[0];
    const double *ARvalue = &ARvalue_[0];

    // Computation
    for (int i = 0; i < ep_count; i++) {
        int iRow = ep_index[i];
        double multi = ep_array[iRow];
        for (int k = ARstart[iRow]; k < ARend_n[iRow]; k++) {
            int index = ARindex[k];
            double value0 = ap_array[index];
            double value1 = value0 + multi * ARvalue[k];
            if (value0 == 0)
                ap_index[ap_count++] = index;
            ap_array[index] = (fabs(value1) < H_TT) ? H_TZ : value1;
        }
    }

    // Try to remove cancellation
    const int apcount1 = ap_count;
    ap_count = 0;
    for (int i = 0; i < apcount1; i++) {
        const int index = ap_index[i];
        const double value = ap_array[index];
        if (fabs(value) > H_TT) {
            ap_index[ap_count++] = index;
        } else {
            ap_array[index] = 0;
        }
    }
    row_ap.count = ap_count;
}

