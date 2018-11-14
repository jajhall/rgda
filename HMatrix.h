#ifndef HMATRIX_H_
#define HMATRIX_H_

#include "HVector.h"

#include <vector>
using namespace std;

class HMatrix {
public:
    void setup(int numCol,
               int numRow,
               const int *Astart,
               const int *Aindex,
               const double *Avalue);
    void scale(double *col_scale, double *row_scale);
    void price_by_col(HVector& row_ap, HVector& row_ep);
    void price_by_row(HVector& row_ap, HVector& row_ep);
    void update(int columnIn, int columnOut);
    double compute_dot(HVector& vector, int iCol);
    void collect_aj(HVector& vector, int iCol, double multi);

    int numCol_;
    int numRow_;
    vector<int> Astart_;
    vector<int> Aindex_;
    vector<double> Avalue_;
private:
    vector<int> ARstart_;
    vector<int> ARend_n_;
    vector<int> ARindex_;
    vector<double> ARvalue_;

};

#endif /* HMATRIX_H_ */
