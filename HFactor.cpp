#include "HFactor.h"
#include "HConst.h"

#include <cmath>
#include <stdexcept>
using namespace std;

void HFactor::setup(const HMatrix *matrix_) {
    // Setup default search parameters
    search_tol = 0.1;
    search_lim = 8;

    // Copy Problem size and (pointer to) coefficient matrix
    numRow = matrix_->numRow_;
    numCol = matrix_->numCol_;
    matrix = matrix_;

    // Allocate for working buffer
    iwork_.assign(numRow + 1, 0);

    // Find Basis matrix limit size
    BlimitX = 0;
    const int *Astart = &matrix->Astart_[0];
    for (int i = 0; i < numCol; i++)
        iwork_[Astart[i + 1] - Astart[i]]++;
    for (int i = numRow, counted = 0; i >= 0 && counted < numRow; i--)
        BlimitX += i * iwork_[i], counted += iwork_[i];
    BlimitX += numRow;

    // Allocate spaces for basis matrix, L, U factor and Update buffer
    int PFlimitp = 1005;
    Bstart_.resize(numRow + 1);
    Lstart_.resize(numRow + 1);
    Ustart_.resize(numRow + 1);
    LRstart_.resize(numRow + 1);
    URstart_.resize(numRow + 1);
    PFstart_.resize(PFlimitp * 2 + 1);
    Bstart_[0] = Lstart_[0] = Ustart_[0] = 0;
    LRstart_[0] = URstart_[0] = PFstart_[0] = 0;

    Bindex_.resize(BlimitX);
    Lindex_.resize(BlimitX * 3);
    Uindex_.resize(BlimitX * 3);
    LRindex_.resize(BlimitX * 3);
    URindex_.resize(BlimitX * 3);
    PFindex_.resize(BlimitX * 4);

    Bvalue_.resize(BlimitX);
    Lvalue_.resize(BlimitX * 3);
    Uvalue_.resize(BlimitX * 3);
    LRvalue_.resize(BlimitX * 3);
    URvalue_.resize(BlimitX * 3);
    PFvalue_.resize(BlimitX * 4);

    // Allocate spaces for pivot records
    permute_.resize(numRow);
    HpivotIndex_.resize(numRow);
    HpivotValue_.resize(numRow);
    PFpivotValue_.resize(PFlimitp);

    // Allocate spaces for Markowitz matrices
    Mcolumn.setup(numRow);

    MCstart.resize(numRow);
    MCcountA.resize(numRow);
    MCplaceN.resize(numRow);
    MCspace.resize(numRow);
    MCminpivot.resize(numRow);

    MRstart.resize(numRow);
    MRcount.resize(numRow);
    MRspace.resize(numRow);
    MRcountb4.resize(numRow);

    MCindex.resize(BlimitX * 2);
    MCvalue.resize(BlimitX * 2);
    MRindex.resize(BlimitX * 2);
    MRfillin.resize(numRow * 8);

    // Allocate spaces for count-link-list
    clinknext.resize(numRow);
    clinklast.resize(numRow);
    clink1.resize(numRow + 1);
    rlinknext.resize(numRow);
    rlinklast.resize(numRow);
    rlink1.resize(numRow + 1);

}

void HFactor::build(int *base) {

    // Prepare build

    build_prepare(base);

    // Build the L, U factor
    build_singleton();
    build_kernel_prepare();
    for (int i = 0; i < nwork; i++) {
        bool good_pivot = build_kernel_search();
        if (good_pivot) {
            build_kernel_eliminate_prepare();
            build_kernel_eliminate_columns();
            build_kernel_eliminate_rows();
        } else {
            throw runtime_error("error-singular-matrix");
        }
    }

    // Finish build
    build_finish(base);

    // Clean up for update
    PFnPivot = 0;
    PFnX = 0;
}

void HFactor::ftran(HVector& vector, double hist_dsty) const {
    const double hyperFTRANL = 0.15;
    const double hyperFTRANU = 0.10;
    const double hyperCANCEL = 0.05;

    // FTRAN Lower
    double curr_dsty = 1.0 * vector.count / numRow;
    if (curr_dsty > hyperCANCEL || hist_dsty > hyperFTRANL) {
        ftran_L(vector);
    } else {
        ftran_LH(vector);
    }

    // FTRAN Update
    ftran_R(vector);
    vector.pack();

    // FTRAN Upper
    curr_dsty = 1.0 * vector.count / numRow;
    if (curr_dsty > hyperCANCEL || hist_dsty > hyperFTRANU) {
        ftran_U(vector);
    } else {
        ftran_UH(vector);
    }
}

void HFactor::btran(HVector& vector, double hist_dsty) const {
    const double hyperBTRANU = 0.15;
    const double hyperBTRANL = 0.10;
    const double hyperCANCEL = 0.05;

    // BTRAN Upper
    double curr_dsty = 1.0 * vector.count / numRow;
    if (curr_dsty > hyperCANCEL || hist_dsty > hyperBTRANU) {
        btran_U(vector);
    } else {
        btran_UH(vector);
    }

    // BTRAN Update
    vector.pack();
    btran_R(vector);

    // BTRAN Lower
    curr_dsty = 1.0 * vector.count / numRow;
    if (curr_dsty > hyperCANCEL || hist_dsty > hyperBTRANL) {
        btran_L(vector);
    } else {
        btran_LH(vector);
    }

}

void HFactor::update(HVector& column, HVector& row_ep, int iRow, bool& re_inv) {
    // 0. See if there's enough spaces for storing the pivots
    int PFlimitp = PFpivotValue_.size();
    if (PFnPivot >= PFlimitp) {
        PFlimitp += 1000;
        PFpivotValue_.resize(PFlimitp);
        PFstart_.resize(PFlimitp * 2 + 1);
    }

    // 1. Translate pivot row to factor index system
    int pivotIndex = permute_[iRow];
    double alpha = column.array[iRow];

    // 2. Check and expand fill space
    const int ucStart = Ustart_[pivotIndex];
    const int ucEnd = Ustart_[pivotIndex + 1];
    const int fillCount = column.packCount + (ucEnd - ucStart + 1)
            + row_ep.packCount;
    int PFMlimitX = PFindex_.size();
    if (PFMlimitX < PFnX + fillCount) {
        PFindex_.resize(PFMlimitX + 4 * BlimitX);
        PFvalue_.resize(PFMlimitX + 4 * BlimitX);
    }

    // 3. Store partial FTRAN result
    copy(&column.packIndex[0], &column.packIndex[0] + column.packCount,
         &PFindex_[PFnX]);
    copy(&column.packValue[0], &column.packValue[0] + column.packCount,
         &PFvalue_[PFnX]);
    PFnX += column.packCount;

    // 4. Store -U_pivotIndex including the pivot entry
    copy(&Uindex_[ucStart], &Uindex_[ucEnd], &PFindex_[PFnX]);
    for (int k = ucStart; k < ucEnd; k++)
        PFvalue_[PFnX++] = -Uvalue_[k];
    PFindex_[PFnX] = iRow;
    PFvalue_[PFnX++] = -HpivotValue_[pivotIndex];
    PFstart_[PFnPivot * 2 + 1] = PFnX;

    // 5. Store partial BTRAN result
    copy(&row_ep.packIndex[0], &row_ep.packIndex[0] + row_ep.packCount,
         &PFindex_[PFnX]);
    copy(&row_ep.packValue[0], &row_ep.packValue[0] + row_ep.packCount,
         &PFvalue_[PFnX]);
    PFnX += row_ep.packCount;
    PFstart_[PFnPivot * 2 + 2] = PFnX;

    // 6. Store pivot element
    PFpivotValue_[PFnPivot++] = alpha;
    int LU_total = LnX + UnX + numRow;
    if (PFnX > LU_total * 3 && PFnPivot >= 50)
        re_inv = true;
}

void HFactor::build_prepare(const int *base) {
    // Reset count
    BnX = LnX = UnX = MCnX = MRnX = npivot = 0;

    // Prepare basis matrix
    const int *Astart = &matrix->Astart_[0];
    const int *Aindex = &matrix->Aindex_[0];
    const double *Avalue = &matrix->Avalue_[0];
    for (int iBColumn = 0; iBColumn < numRow; iBColumn++) {
        int iAColumn = base[iBColumn];
        if (iAColumn >= numCol) { // Logical column
            Bindex_[BnX] = iAColumn - numCol;
            Bvalue_[BnX++] = 1;
        } else { // Structural column
            const int start = Astart[iAColumn];
            const int end = Astart[iAColumn + 1];
            for (int k = start; k < end; k++) {
                Bindex_[BnX] = Aindex[k];
                Bvalue_[BnX++] = Avalue[k];
            }
        }
        Bstart_[iBColumn + 1] = BnX;
    }

    // Prepare available column index
    nwork = numRow;
    for (int i = 0; i < numRow; i++)
        iwork_[i] = i;

    // Prepare basis row count
    MRcountb4.assign(numRow, 0);
    for (int k = 0; k < BnX; k++)
        MRcountb4[Bindex_[k]]++;
}

void HFactor::build_singleton() {
    // Special singletons: all column with 1 will be consider as slacks
    for (int iAvail = 0; iAvail < nwork; iAvail++) {
        pivotCol = iwork_[iAvail];
        const int start = Bstart_[pivotCol];
        const int count = Bstart_[pivotCol + 1] - start;
        if (count == 1 && Bvalue_[start] == 1.0) {
            pivotRow = Bindex_[start];
            iwork_[iAvail] = -1;
            MRcountb4[pivotRow] = 0;
            permute_[pivotCol] = pivotRow;
            Lstart_[npivot + 1] = LnX;
            Ustart_[npivot + 1] = UnX;
            HpivotIndex_[npivot] = pivotRow;
            HpivotValue_[npivot++] = 1;
        }
    }

    // Clear up slack part
    if (npivot) {
        int colSavedCount = 0;
        for (int iAvail = 0; iAvail < nwork; iAvail++)
            if (iwork_[iAvail] != -1)
                iwork_[colSavedCount++] = iwork_[iAvail];
        nwork = colSavedCount;
    }

    // Now move on to real singletons
    while (nwork > 0) {
        for (int iAvail = 0; iAvail < nwork; iAvail++) {
            pivotCol = iwork_[iAvail];
            const int start = Bstart_[pivotCol];
            const int end = Bstart_[pivotCol + 1];
            int found_row_singleton = 0;
            int columnCount = 0;

            // Search for singleton
            for (int k = start; k < end; k++) {
                const int iRow = Bindex_[k];
                if (MRcountb4[iRow] == 1) {
                    pivotK = k;
                    found_row_singleton = 1;
                    break;
                }
                if (MRcountb4[iRow] > 1) {
                    pivotK = k;
                    columnCount++;
                }
            }

            // Record the singleton
            if (found_row_singleton || columnCount == 1) {
                // Disable this column
                iwork_[iAvail] = -1;
                // Record to L, U factor
                if (found_row_singleton)
                    build_singleton_row();
                else
                    build_singleton_column();
            }
        }

        // Shrink the available column list
        int colSavedCount = 0;
        for (int iAvail = 0; iAvail < nwork; iAvail++)
            if (iwork_[iAvail] != -1)
                iwork_[colSavedCount++] = iwork_[iAvail];

        // No singleton found in the last pass
        if (colSavedCount == nwork)
            break;

        // Move to the next search loop
        nwork = colSavedCount;
    }

}

void HFactor::build_singleton_column() {
    // Setup pointer
    const int start = Bstart_[pivotCol];
    const int end = Bstart_[pivotCol + 1];

    // Fill first part to U
    for (int k = start; k < pivotK; k++) {
        Uindex_[UnX] = Bindex_[k];
        Uvalue_[UnX++] = Bvalue_[k];
    }

    // Fill second part to U
    for (int k = pivotK + 1; k < end; k++) {
        Uindex_[UnX] = Bindex_[k];
        Uvalue_[UnX++] = Bvalue_[k];
    }

    // Record pivot
    pivotRow = Bindex_[pivotK];
    MRcountb4[pivotRow] = 0;
    permute_[pivotCol] = pivotRow;
    Lstart_[npivot + 1] = LnX;
    Ustart_[npivot + 1] = UnX;
    HpivotIndex_[npivot] = pivotRow;
    HpivotValue_[npivot++] = Bvalue_[pivotK];
}

void HFactor::build_singleton_row() {
    // Pivot value inverse
    const double pivotX_inv = 1 / Bvalue_[pivotK];

    // Setup pointer
    const int start = Bstart_[pivotCol];
    const int end = Bstart_[pivotCol + 1];

    // Fill first part
    for (int k = start; k < pivotK; k++) {
        const int iRow = Bindex_[k];
        if (MRcountb4[iRow] > 0) {
            Lindex_[LnX] = iRow;
            Lvalue_[LnX++] = Bvalue_[k] * pivotX_inv;
        } else {
            Uindex_[UnX] = iRow;
            Uvalue_[UnX++] = Bvalue_[k];
        }
        MRcountb4[iRow]--;
    }

    // Fill second part
    for (int k = pivotK + 1; k < end; k++) {
        const int iRow = Bindex_[k];
        if (MRcountb4[iRow] > 0) {
            Lindex_[LnX] = iRow;
            Lvalue_[LnX++] = Bvalue_[k] * pivotX_inv;
        } else {
            Uindex_[UnX] = iRow;
            Uvalue_[UnX++] = Bvalue_[k];
        }
        MRcountb4[iRow]--;
    }

    // Record pivot
    pivotRow = Bindex_[pivotK];
    MRcountb4[pivotRow] = 0;
    permute_[pivotCol] = pivotRow;
    Lstart_[npivot + 1] = LnX;
    Ustart_[npivot + 1] = UnX;
    HpivotIndex_[npivot] = pivotRow;
    HpivotValue_[npivot++] = Bvalue_[pivotK];
}

void HFactor::build_kernel_prepare() {
    // Reset link list
    clink1.assign(numRow + 1, -1);
    rlink1.assign(numRow + 1, -1);

    // Prepare kernel row pointer, row link list
    for (int iRow = 0; iRow < numRow; iRow++) {
        const int count = MRcountb4[iRow];
        if (count > 0) {
            // Kernel row pointer
            MRstart[iRow] = MRnX;
            MRcount[iRow] = 0;
            MRspace[iRow] = count * 2;
            MRnX += MRspace[iRow];
            rlinkAdd(iRow, count);
        }
    }

    // Prepare column pointer, kernel matrix, column link list
    for (int iAvail = 0; iAvail < nwork; iAvail++) {
        const int iCol = iwork_[iAvail];
        const int length = Bstart_[iCol + 1] - Bstart_[iCol];

        // Kernel column pointer
        MCstart[iCol] = MCnX;
        MCcountA[iCol] = 0;
        MCplaceN[iCol] = length * 2;
        MCspace[iCol] = length * 2;
        MCnX += MCspace[iCol];

        // Kernel matrix element
        double maxValue = 0;
        for (int k = Bstart_[iCol]; k < Bstart_[iCol + 1]; k++) {
            const int iRow = Bindex_[k];
            const double value = Bvalue_[k];
            if (MRcountb4[iRow] > 0) {
                // Active element
                int row_put = MRstart[iRow] + MRcount[iRow]++;
                int col_put = MCstart[iCol] + MCcountA[iCol]++;
                MRindex[row_put] = iCol;
                MCindex[col_put] = iRow;
                MCvalue[col_put] = value;
                maxValue = max(maxValue, fabs(value));
            } else {
                // Non active element
                int non_put = MCstart[iCol] + (--MCplaceN[iCol]);
                MCindex[non_put] = iRow;
                MCvalue[non_put] = value;
            }
        }
        MCminpivot[iCol] = maxValue * search_tol;
        // Add column index to link list
        clinkAdd(iCol, MCcountA[iCol]);
    }
}

bool HFactor::build_kernel_search() {
    int searchCount = 0;
    const float M_limit = 1.0 * numRow * numRow;
    float M_best = M_limit;

    // Pick a column singleton
    if (clink1[1] != -1) {
        pivotCol = clink1[1];
        pivotRow = MCindex[MCstart[pivotCol]];
        pivotK = MCstart[pivotCol];
        return true;
    }

    // Pick a row singleton
    if (rlink1[1] != -1) {
        pivotRow = rlink1[1];
        pivotCol = MRindex[MRstart[pivotRow]];
        pivotK = MCstart[pivotCol];
        while (MCindex[pivotK] != pivotRow)
            pivotK++;
        return true;
    }

    for (int count = 2; count <= numRow; count++) {
        // Check column lists
        for (int iCol = clink1[count]; iCol != -1; iCol = clinknext[iCol]) {
            const double minpivot = MCminpivot[iCol];
            const int start = MCstart[iCol], end = start + MCcountA[iCol];
            for (int k = start; k < end; k++) {
                if (fabs(MCvalue[k]) >= minpivot) {
                    const int iRow = MCindex[k];
                    const int rowCount = MRcount[iRow];
                    const float M_value = 1.0 * (count - 1) * (rowCount - 1);
                    if (M_value < M_best) {
                        M_best = M_value;
                        pivotCol = iCol;
                        pivotRow = iRow;
                        pivotK = k;
                        if (rowCount < count) // MM <= (count-1)^2
                            return true;
                    }
                }
            }
            if (searchCount++ >= search_lim && M_best < M_limit)
                return true;
        }

        // Check row lists
        for (int iRow = rlink1[count]; iRow != -1; iRow = rlinknext[iRow]) {
            const int start = MRstart[iRow], end = start + MRcount[iRow];
            for (int k = start; k < end; k++) {
                int iCol = MRindex[k];
                int columnCount = MCcountA[iCol];
                float M_value = 1.0 * (count - 1) * (columnCount - 1);
                if (M_value < M_best) {
                    int ifind = MCstart[iCol];
                    while (MCindex[ifind] != iRow)
                        ifind++;
                    if (fabs(MCvalue[ifind]) >= MCminpivot[iCol]) {
                        M_best = M_value;
                        pivotCol = iCol;
                        pivotRow = iRow;
                        pivotK = ifind;
                        if (columnCount <= count) // MM <= count * (count-1)
                            return true;
                    }
                }
            }
            if (searchCount++ >= search_lim && M_best < M_limit)
                return true;
        }

    }

    // Need to check if is singular
    return (M_best < M_limit);
}

void HFactor::build_kernel_eliminate_prepare() {
    // Prepare pivot column
    const double pivotX_inv = MCvalue[pivotK];
    const int startA = MCstart[pivotCol];
    const int endA = startA + MCcountA[pivotCol];

    int McolumnCount = 0;
    int *McolumnIndex = &Mcolumn.index[0];
    char *McolumnMark = (char *) &Mcolumn.iwork[0];
    double *McolumnArray = &Mcolumn.array[0];
    for (int k = startA; k < endA; k++) {
        if (pivotK != k) {
            const int iRow = MCindex[k];
            const double value = MCvalue[k] / pivotX_inv;
            // Prepare pivot column elimination buffer
            McolumnIndex[McolumnCount++] = iRow;
            McolumnArray[iRow] = value;
            McolumnMark[iRow] = 1;

            // Fill factor L
            Lindex_[LnX] = iRow;
            Lvalue_[LnX++] = value;
            // Remove from row matrix
            MRcountb4[iRow] = MRcount[iRow];
            int iRemove = MRstart[iRow];
            int iTake = iRemove + (--MRcount[iRow]);
            while (MRindex[iRemove] != pivotCol)
                iRemove++;
            MRindex[iRemove] = MRindex[iTake];
        }
    }
    Mcolumn.count = McolumnCount;

    // Fill factor U
    const int startN = startA + MCplaceN[pivotCol];
    const int endN = startA + MCspace[pivotCol];
    copy(&MCindex[startN], &MCindex[endN], &Uindex_[UnX]);
    copy(&MCvalue[startN], &MCvalue[endN], &Uvalue_[UnX]);
    UnX += endN - startN;

    // Record pivot
    permute_[pivotCol] = pivotRow;
    Lstart_[npivot + 1] = LnX;
    Ustart_[npivot + 1] = UnX;
    HpivotIndex_[npivot] = pivotRow;
    HpivotValue_[npivot++] = MCvalue[pivotK];
    clinkDel(pivotCol);
    rlinkDel(pivotRow);

    // Check factor L, U space for future
    const int LlimitX = Lindex_.size();
    if (LlimitX < LnX + numRow) {
        Lindex_.resize(LlimitX + BlimitX * 3);
        Lvalue_.resize(LlimitX + BlimitX * 3);
    }
    const int UlimitX = Uindex_.size();
    if (UlimitX < UnX + numRow) {
        Uindex_.resize(UlimitX + BlimitX * 3);
        Uvalue_.resize(UlimitX + BlimitX * 3);
    }
}

void HFactor::build_kernel_eliminate_columns() {
    // Alias to M column
    const int McolumnCount = Mcolumn.count;
    const int *McolumnIndex = &Mcolumn.index[0];
    char *McolumnMark = (char *) &Mcolumn.iwork[0];
    const double *McolumnArray = &Mcolumn.array[0];

    // Loop over pivot row to eliminate other column
    MRnFillin = 0;
    const int PRstart = MRstart[pivotRow];
    const int PRend = PRstart + MRcount[pivotRow];
    for (int iPR = PRstart; iPR < PRend; iPR++) {
        int iCol = MRindex[iPR];
        if (iCol == pivotCol)
            continue;

        // Setup updated column pointer
        const int countA_before = MCcountA[iCol];
        const int startA = MCstart[iCol];
        const int endA = startA + (--MCcountA[iCol]);
        const int startN = startA + (--MCplaceN[iCol]);

        // Move pivot to non active part (NOTE sequence!!)
        int ipivot = startA;
        while (MCindex[ipivot] != pivotRow)
            ipivot++;
        const double pivotX = MCvalue[ipivot];
        MCindex[ipivot] = MCindex[endA];
        MCvalue[ipivot] = MCvalue[endA];
        MCindex[startN] = pivotRow;
        MCvalue[startN] = pivotX;

        // Elimination on the overlapping part
        int count_fillin = McolumnCount;
        int count_cancel = 0;
        double maxValue = 0;
        for (int k = startA; k < endA; k++) {
            int iRow = MCindex[k];
            double value = MCvalue[k];
            if (McolumnMark[iRow]) {
                McolumnMark[iRow] = 0;
                count_fillin--;

                value -= pivotX * McolumnArray[iRow];
                if (fabs(value) < H_TT) {
                    value = 0;
                    count_cancel++;
                }
                MCvalue[k] = value;
            }
            maxValue = max(maxValue, fabs(value));
        }

        // Remove cancellation gaps
        if (count_cancel > 0) {
            int countA = startA;
            for (int i = startA; i < endA; i++) {
                if (MCvalue[i] != 0) {
                    // Shrink the cancellation gap in the column matrix
                    MCindex[countA] = MCindex[i];
                    MCvalue[countA++] = MCvalue[i];
                } else {
                    // Remove cancellation from the row matrix
                    int iRow = MCindex[i];
                    int iRemove = MRstart[iRow];
                    int iTake = iRemove + (--MRcount[iRow]);
                    while (MRindex[iRemove] != iCol)
                        iRemove++;
                    MRindex[iRemove] = MRindex[iTake];
                }
            }
            MCcountA[iCol] = countA - startA;
        }

        // Prepare space for fill in
        if (MCplaceN[iCol] - MCcountA[iCol] < count_fillin)
            build_expand_column(iCol, count_fillin);

        // Fill into current column
        if (count_fillin > 0) {
            int ifillin = MCstart[iCol] + MCcountA[iCol];
            for (int i = 0; i < McolumnCount; i++) {
                int iRow = McolumnIndex[i];
                if (McolumnMark[iRow]) {
                    const double value = -pivotX * McolumnArray[iRow];
                    maxValue = max(maxValue, fabs(value));
                    MCindex[ifillin] = iRow;
                    MCvalue[ifillin++] = value;
                    MRfillin[MRnFillin++] = iRow;
                    MRfillin[MRnFillin++] = iCol;
                } else
                    McolumnMark[iRow] = 1;
            }
            MCcountA[iCol] = ifillin - MCstart[iCol];

            // Check MR fill-in record space for future
            const int MRfillinLimit = MRfillin.size();
            if (MRfillinLimit < MRnFillin + numRow * 2) {
                MRfillin.resize(MRfillinLimit + numRow * 4);
            }

        } else {
            // Restore mark when there's no fill-in
            for (int i = 0; i < McolumnCount; i++)
                McolumnMark[McolumnIndex[i]] = 1;
        }

        // Fix max value and link list
        MCminpivot[iCol] = maxValue * search_tol;
        if (countA_before != MCcountA[iCol]) {
            clinkDel(iCol);
            clinkAdd(iCol, MCcountA[iCol]);
        }
    }

    // Restore pivot column mark
    for (int i = 0; i < McolumnCount; i++)
        McolumnMark[McolumnIndex[i]] = 0;

}

void HFactor::build_kernel_eliminate_rows() {
    // Fill in to row matrix
    for (int i = 0; i < MRnFillin; i += 2) {
        const int iRow = MRfillin[i];
        const int iCol = MRfillin[i + 1];
        if (MRcount[iRow] == MRspace[iRow])
            build_expand_row(iRow);
        const int ifillin = MRstart[iRow] + MRcount[iRow]++;
        MRindex[ifillin] = iCol;
    }

    // Update link list for the all of the pivot row index
    const int PCstartA = MCstart[pivotCol];
    const int PCendA = PCstartA + MCcountA[pivotCol];
    for (int i = PCstartA; i < PCendA; i++) {
        int iRow = MCindex[i];
        if (iRow != pivotRow) {
            if (MRcountb4[iRow] != MRcount[iRow]) {
                rlinkDel(iRow);
                rlinkAdd(iRow, MRcount[iRow]);
            }
        }
    }
}

void HFactor::build_expand_column(const int iCol, const int fillin) {
    // Setup new space
    const int space_inc = max(MCspace[iCol], fillin);
    const int space_new = MCspace[iCol] + space_inc;
    // Expand the column matrix if needed
    const int MClimitX = MCindex.size();
    if (MClimitX <= space_new + MCnX) {
        MCindex.resize(MClimitX + BlimitX * 2);
        MCvalue.resize(MClimitX + BlimitX * 2);
    }
    // Copy column while expand the gap
    const int startA = MCstart[iCol];
    const int startN = startA + MCplaceN[iCol];
    const int endA = startA + MCcountA[iCol];
    const int endN = startA + MCspace[iCol];
    const int placeN_new = space_inc + MCplaceN[iCol];
    const int startN_new = MCnX + placeN_new;
    copy(&MCindex[startA], &MCindex[endA], &MCindex[MCnX]);
    copy(&MCvalue[startA], &MCvalue[endA], &MCvalue[MCnX]);
    copy(&MCindex[startN], &MCindex[endN], &MCindex[startN_new]);
    copy(&MCvalue[startN], &MCvalue[endN], &MCvalue[startN_new]);
    // Record the new pointers
    MCstart[iCol] = MCnX;
    MCplaceN[iCol] = placeN_new;
    MCspace[iCol] = space_new;
    MCnX += space_new;
}

void HFactor::build_expand_row(const int iRow) {
    // Setup new space
    const int space_new = MRspace[iRow] * 2;
    // Expand row matrix if needed
    const int MRlimitX = MRindex.size();
    if (MRlimitX <= space_new + MRnX) {
        MRindex.resize(MRlimitX + BlimitX);
    }
    // Copy row matrix
    const int start = MRstart[iRow];
    const int end = start + MRcount[iRow];
    copy(&MRindex[start], &MRindex[end], &MRindex[MRnX]);
    // Setup new pointer
    MRstart[iRow] = MRnX;
    MRspace[iRow] = space_new;
    MRnX += space_new;
}

void HFactor::build_finish(int *base) {
    // Finish: base index
    iwork_.assign(base, base + numRow);
    for (int i = 0; i < numRow; i++)
        base[permute_[i]] = iwork_[i];

    for (int i = 0; i < numRow; i++)
        permute_[HpivotIndex_[i]] = i;

    // Finish: LR factor
    MRcountb4.assign(numRow, 0);
    LRstart_[0] = 0;
    if (LRindex_.size() != Lindex_.size()) {
        LRindex_.resize(Lindex_.size());
        LRvalue_.resize(Lvalue_.size());
    }

    for (int k = 0; k < LnX; k++)
        MRcountb4[permute_[Lindex_[k]]]++;
    for (int i = 1; i <= numRow; i++)
        LRstart_[i] = LRstart_[i - 1] + MRcountb4[i - 1];
    MRcountb4.assign(LRstart_.begin(), LRstart_.begin() + numRow);

    for (int i = 0; i < numRow; i++) {
        int iCol = HpivotIndex_[i];
        for (int k = Lstart_[i]; k < Lstart_[i + 1]; k++) {
            int iRow = permute_[Lindex_[k]];
            int iPut = MRcountb4[iRow]++;
            LRindex_[iPut] = iCol;
            LRvalue_[iPut] = Lvalue_[k];
        }
    }

    // Finish: UR factor
    MRcountb4.assign(numRow, 0);
    URstart_[0] = 0;
    if (URindex_.size() != Uindex_.size()) {
        URindex_.resize(Uindex_.size());
        URvalue_.resize(Uindex_.size());
    }

    for (int k = 0; k < UnX; k++)
        MRcountb4[permute_[Uindex_[k]]]++;
    for (int i = 1; i <= numRow; i++)
        URstart_[i] = URstart_[i - 1] + MRcountb4[i - 1];
    MRcountb4.assign(URstart_.begin(), URstart_.begin() + numRow);

    for (int i = 0; i < numRow; i++) {
        int iCol = HpivotIndex_[i];
        for (int k = Ustart_[i]; k < Ustart_[i + 1]; k++) {
            int iRow = permute_[Uindex_[k]];
            int iPut = MRcountb4[iRow]++;
            URindex_[iPut] = iCol;
            URvalue_[iPut] = Uvalue_[k];
        }
    }
}

void HFactor::ftran_L(HVector& vector) const {
    const int Lsize = numRow;
    const int *Lstart = &Lstart_[0];
    const int *Lindex = &Lindex_[0];
    const double *Lvalue = &Lvalue_[0];
    const int *LpivotIndex = &HpivotIndex_[0];

    int vector_count = 0;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    for (int i = 0; i < Lsize; i++) {
        int pivotRow = LpivotIndex[i];
        const double pivotX = vector_array[pivotRow];
        if (fabs(pivotX) > H_TT) {
            vector_index[vector_count++] = pivotRow;
            for (int k = Lstart[i]; k < Lstart[i + 1]; k++)
                vector_array[Lindex[k]] -= pivotX * Lvalue[k];
        } else
            vector_array[pivotRow] = 0;
    }
    vector.count = vector_count;
}

void HFactor::ftran_R(HVector& vector) const {
    const int *PFstart = &PFstart_[0];
    const int *PFindex = &PFindex_[0];
    const double *PFvalue = &PFvalue_[0];
    const double *PFpivotValue = &PFpivotValue_[0];

    int vector_count = vector.count;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    for (int i = 0; i < PFnPivot; i++) {
        const int Xstart = PFstart[i * 2];
        const int Xend = PFstart[i * 2 + 1];
        const int Ystart = Xend;
        const int Yend = PFstart[i * 2 + 2];

        double pivotX = 0;
        for (int k = Ystart; k < Yend; k++)
            pivotX += PFvalue[k] * vector_array[PFindex[k]];

        if (fabs(pivotX) > H_TT) {
            pivotX /= PFpivotValue[i];
            for (int k = Xstart; k < Xend; k++) {
                const int index = PFindex[k];
                const double value0 = vector_array[index];
                const double value1 = value0 - pivotX * PFvalue[k];
                if (value0 == 0)
                    vector_index[vector_count++] = index;
                vector_array[index] = (fabs(value1) < H_TT) ? H_TZ : value1;
            }
        }
    }

    // Remove cancellation
    const int vector_count1 = vector_count;
    vector_count = 0;
    for (int i = 0; i < vector_count1; i++) {
        const int index = vector_index[i];
        const double value = vector_array[index];
        if (fabs(value) > H_TT) {
            vector_index[vector_count++] = index;
        } else {
            vector_array[index] = 0;
        }
    }
    vector.count = vector_count;

}

void HFactor::ftran_U(HVector& vector) const {
    const int Usize = numRow;
    const int *Ustart = &Ustart_[0];
    const int *Uindex = &Uindex_[0];
    const double *Uvalue = &Uvalue_[0];
    const int *UpivotIndex = &HpivotIndex_[0];
    const double *UpivotValue = &HpivotValue_[0];

    int vector_count = 0;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    for (int i = Usize - 1; i >= 0; i--) {
        int pivotRow = UpivotIndex[i];
        double pivotX = vector_array[pivotRow];
        if (fabs(pivotX) > H_TT) {
            pivotX /= UpivotValue[i];
            vector_index[vector_count++] = pivotRow;
            vector_array[pivotRow] = pivotX;
            for (int k = Ustart[i]; k < Ustart[i + 1]; k++)
                vector_array[Uindex[k]] -= pivotX * Uvalue[k];
        } else
            vector_array[pivotRow] = 0;
    }
    vector.count = vector_count;
}

void HFactor::btran_U(HVector& vector) const {
    const int URsize = numRow;
    const int *URstart = &URstart_[0];
    const int *URindex = &URindex_[0];
    const double *URvalue = &URvalue_[0];
    const int *URpivotIndex = &HpivotIndex_[0];
    const double *URpivotValue = &HpivotValue_[0];

    int vector_count = 0;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    for (int i = 0; i < URsize; i++) {
        int pivotRow = URpivotIndex[i];
        double pivotX = vector_array[pivotRow];
        if (fabs(pivotX) > H_TT) {
            pivotX /= URpivotValue[i];
            vector_index[vector_count++] = pivotRow;
            vector_array[pivotRow] = pivotX;
            for (int k = URstart[i]; k < URstart[i + 1]; k++)
                vector_array[URindex[k]] -= pivotX * URvalue[k];
        } else
            vector_array[pivotRow] = 0;
    }
    vector.count = vector_count;
}

void HFactor::btran_R(HVector& vector) const {
    const int *PFstart = &PFstart_[0];
    const int *PFindex = &PFindex_[0];
    const double *PFvalue = &PFvalue_[0];
    const double *PFpivotValue = &PFpivotValue_[0];

    int vector_count = vector.count;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    for (int i = PFnPivot - 1; i >= 0; i--) {
        const int Xstart = PFstart[i * 2];
        const int Xend = PFstart[i * 2 + 1];
        const int Ystart = Xend;
        const int Yend = PFstart[i * 2 + 2];

        double pivotX = 0;
        for (int k = Xstart; k < Xend; k++)
            pivotX += PFvalue[k] * vector_array[PFindex[k]];

        if (fabs(pivotX) > H_TT) {
            pivotX /= PFpivotValue[i];
            for (int k = Ystart; k < Yend; k++) {
                const int index = PFindex[k];
                const double value0 = vector_array[index];
                const double value1 = value0 - pivotX * PFvalue[k];
                if (value0 == 0)
                    vector_index[vector_count++] = index;
                vector_array[index] = (fabs(value1) < H_TT) ? H_TZ : value1;
            }
        }
    }

    // Remove cancellation
    const int vector_count1 = vector_count;
    vector_count = 0;
    for (int i = 0; i < vector_count1; i++) {
        const int index = vector_index[i];
        const double value = vector_array[index];
        if (fabs(value) > H_TT) {
            vector_index[vector_count++] = index;
        } else {
            vector_array[index] = 0;
        }
    }
    vector.count = vector_count;
}

void HFactor::btran_L(HVector& vector) const {
    const int LRsize = numRow;
    const int *LRstart = &LRstart_[0];
    const int *LRindex = &LRindex_[0];
    const double *LRvalue = &LRvalue_[0];
    const int *LRpivotIndex = &HpivotIndex_[0];

    int vector_count = 0;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    for (int i = LRsize - 1; i >= 0; i--) {
        int pivotRow = LRpivotIndex[i];
        const double pivotX = vector_array[pivotRow];
        if (fabs(pivotX) > H_TT) {
            vector_index[vector_count++] = pivotRow;
            vector_array[pivotRow] = pivotX;
            for (int k = LRstart[i]; k < LRstart[i + 1]; k++)
                vector_array[LRindex[k]] -= pivotX * LRvalue[k];
        } else
            vector_array[pivotRow] = 0;
    }
    vector.count = vector_count;
}

void HFactor::ftran_LH(HVector& vector) const {
    // Alias to L
    const int Lsize = numRow;
    const int *Lstart = &Lstart_[0];
    const int *Lindex = &Lindex_[0];
    const double *Lvalue = &Lvalue_[0];
    const int *Lpermute = &permute_[0];
    const int *LpivotIndex = &HpivotIndex_[0];

    int vector_count = vector.count;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    // 1. Build list
    const int *Lend = &Lstart[1];
    char *listMark = (char *) &vector.iwork[0];
    int *listIndex = &vector.iwork[Lsize];
    int *listStack = &vector.iwork[Lsize * 2];
    int listCount = 0;

    for (int i = 0; i < vector_count; i++) {
        // Skip touched index
        int iTrans = Lpermute[vector_index[i]];
        if (listMark[iTrans])
            continue;

        int Hi = iTrans; // H matrix pivot index
        int Hk = Lstart[Hi]; // H matrix non zero position
        int nStack = -1; // Usage of the stack (-1 not used)

        listMark[Hi] = 1; // Mark this as touched

        for (;;) {
            if (Hk < Lend[Hi]) {
                int Hi_sub = Lpermute[Lindex[Hk++]];
                if (listMark[Hi_sub] == 0) { // Go to a child
                    listMark[Hi_sub] = 1; // Mark as touched
                    listStack[++nStack] = Hi; // Store current into stack
                    listStack[++nStack] = Hk;
                    Hi = Hi_sub; // Replace current with child
                    Hk = Lstart[Hi];
                }
            } else {
                listIndex[listCount++] = Hi;
                if (nStack == -1) // Quit on empty stack
                    break;
                Hk = listStack[nStack--]; // Back to last in stack
                Hi = listStack[nStack--];
            }
        }
    }

    // 2. Solve with list
    vector_count = 0;
    for (int iList = listCount - 1; iList >= 0; iList--) {
        int i = listIndex[iList];
        listMark[i] = 0;
        int pivotRow = LpivotIndex[i];
        double pivotX = vector_array[pivotRow];
        if (fabs(pivotX) > H_TT) {
            vector_index[vector_count++] = pivotRow;
            for (int k = Lstart[i]; k < Lstart[i + 1]; k++)
                vector_array[Lindex[k]] -= pivotX * Lvalue[k];
        } else
            vector_array[pivotRow] = 0;
    }
    vector.count = vector_count;
}

void HFactor::ftran_UH(HVector& vector) const {
    // Alias to U
    const int Usize = numRow;
    const int *Ustart = &Ustart_[0];
    const int *Uindex = &Uindex_[0];
    const double *Uvalue = &Uvalue_[0];
    const int *Upermute = &permute_[0];
    const int *UpivotIndex = &HpivotIndex_[0];
    const double *UpivotValue = &HpivotValue_[0];

    int vector_count = vector.count;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    // 1. Build list
    const int *Lend = &Ustart[1];
    char *listMark = (char *) &vector.iwork[0];
    int *listIndex = &vector.iwork[Usize];
    int *listStack = &vector.iwork[Usize * 2];
    int listCount = 0;

    for (int i = 0; i < vector_count; i++) {
        // Skip touched index
        int iTrans = Upermute[vector_index[i]];
        if (listMark[iTrans])
            continue;

        int Hi = iTrans; // H matrix pivot index
        int Hk = Ustart[Hi]; // H matrix non zero position
        int nStack = -1; // Usage of the stack (-1 not used)

        listMark[Hi] = 1; // Mark this as touched

        for (;;) {
            if (Hk < Lend[Hi]) {
                int Hi_sub = Upermute[Uindex[Hk++]];
                if (listMark[Hi_sub] == 0) { // Go to a child
                    listMark[Hi_sub] = 1; // Mark as touched
                    listStack[++nStack] = Hi; // Store current into stack
                    listStack[++nStack] = Hk;
                    Hi = Hi_sub; // Replace current with child
                    Hk = Ustart[Hi];
                }
            } else {
                listIndex[listCount++] = Hi;
                if (nStack == -1) // Quit on empty stack
                    break;
                Hk = listStack[nStack--]; // Back to last in stack
                Hi = listStack[nStack--];
            }
        }
    }

    // 2. Solve with list
    vector_count = 0;
    for (int iList = listCount - 1; iList >= 0; iList--) {
        int i = listIndex[iList];
        listMark[i] = 0;
        int pivotRow = UpivotIndex[i];
        double pivotX = vector_array[pivotRow];
        if (fabs(pivotX) > H_TT) {
            pivotX /= UpivotValue[i];
            vector_index[vector_count++] = pivotRow;
            vector_array[pivotRow] = pivotX;
            for (int k = Ustart[i]; k < Ustart[i + 1]; k++)
                vector_array[Uindex[k]] -= pivotX * Uvalue[k];
        } else
            vector_array[pivotRow] = 0;
    }
    vector.count = vector_count;
}

void HFactor::btran_UH(HVector& vector) const {
    // Alias to U
    const int URsize = numRow;
    const int *URstart = &URstart_[0];
    const int *URindex = &URindex_[0];
    const double *URvalue = &URvalue_[0];
    const int *URpermute = &permute_[0];
    const int *URpivotIndex = &HpivotIndex_[0];
    const double *URpivotValue = &HpivotValue_[0];

    int vector_count = vector.count;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    // 1. Build list
    const int *Lend = &URstart[1];
    char *listMark = (char *) &vector.iwork[0];
    int *listIndex = &vector.iwork[URsize];
    int *listStack = &vector.iwork[URsize * 2];
    int listCount = 0;

    for (int i = 0; i < vector_count; i++) {
        // Skip touched index
        int iTrans = URpermute[vector_index[i]];
        if (listMark[iTrans])
            continue;

        int Hi = iTrans; // H matrix pivot index
        int Hk = URstart[Hi]; // H matrix non zero position
        int nStack = -1; // Usage of the stack (-1 not used)

        listMark[Hi] = 1; // Mark this as touched

        for (;;) {
            if (Hk < Lend[Hi]) {
                int Hi_sub = URpermute[URindex[Hk++]];
                if (listMark[Hi_sub] == 0) { // Go to a child
                    listMark[Hi_sub] = 1; // Mark as touched
                    listStack[++nStack] = Hi; // Store current into stack
                    listStack[++nStack] = Hk;
                    Hi = Hi_sub; // Replace current with child
                    Hk = URstart[Hi];
                }
            } else {
                listIndex[listCount++] = Hi;
                if (nStack == -1) // Quit on empty stack
                    break;
                Hk = listStack[nStack--]; // Back to last in stack
                Hi = listStack[nStack--];
            }
        }
    }

    // 2. Solve with list
    vector_count = 0;
    for (int iList = listCount - 1; iList >= 0; iList--) {
        int i = listIndex[iList];
        listMark[i] = 0;
        int pivotRow = URpivotIndex[i];
        double pivotX = vector_array[pivotRow];
        if (fabs(pivotX) > H_TT) {
            pivotX /= URpivotValue[i];
            vector_index[vector_count++] = pivotRow;
            vector_array[pivotRow] = pivotX;
            for (int k = URstart[i]; k < URstart[i + 1]; k++)
                vector_array[URindex[k]] -= pivotX * URvalue[k];
        } else
            vector_array[pivotRow] = 0;
    }
    vector.count = vector_count;
}

void HFactor::btran_LH(HVector& vector) const {
    // Alias to LR
    const int LRsize = numRow;
    const int *LRstart = &LRstart_[0];
    const int *LRindex = &LRindex_[0];
    const double *LRvalue = &LRvalue_[0];
    const int *LRpermute = &permute_[0];
    const int *LRpivotIndex = &HpivotIndex_[0];

    int vector_count = vector.count;
    int *vector_index = &vector.index[0];
    double *vector_array = &vector.array[0];

    // 1. Build list
    const int *Lend = &LRstart[1];
    char *listMark = (char *) &vector.iwork[0];
    int *listIndex = &vector.iwork[LRsize];
    int *listStack = &vector.iwork[LRsize * 2];
    int listCount = 0;

    for (int i = 0; i < vector_count; i++) {
        // Skip touched index
        int iTrans = LRpermute[vector_index[i]];
        if (listMark[iTrans])
            continue;

        int Hi = iTrans; // H matrix pivot index
        int Hk = LRstart[Hi]; // H matrix non zero position
        int nStack = -1; // Usage of the stack (-1 not used)

        listMark[Hi] = 1; // Mark this as touched

        for (;;) {
            if (Hk < Lend[Hi]) {
                int Hi_sub = LRpermute[LRindex[Hk++]];
                if (listMark[Hi_sub] == 0) { // Go to a child
                    listMark[Hi_sub] = 1; // Mark as touched
                    listStack[++nStack] = Hi; // Store current into stack
                    listStack[++nStack] = Hk;
                    Hi = Hi_sub; // Replace current with child
                    Hk = LRstart[Hi];
                }
            } else {
                listIndex[listCount++] = Hi;
                if (nStack == -1) // Quit on empty stack
                    break;
                Hk = listStack[nStack--]; // Back to last in stack
                Hi = listStack[nStack--];
            }
        }
    }

    // 2. Solve with list
    vector_count = 0;
    for (int iList = listCount - 1; iList >= 0; iList--) {
        int i = listIndex[iList];
        listMark[i] = 0;
        int pivotRow = LRpivotIndex[i];
        double pivotX = vector_array[pivotRow];
        if (fabs(pivotX) > H_TT) {
            vector_index[vector_count++] = pivotRow;
            for (int k = LRstart[i]; k < LRstart[i + 1]; k++)
                vector_array[LRindex[k]] -= pivotX * LRvalue[k];
        } else
            vector_array[pivotRow] = 0;
    }
    vector.count = vector_count;
}
