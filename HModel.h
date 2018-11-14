#ifndef HMODEL_H_
#define HMODEL_H_

#include "HInput.h"
#include "HVector.h"
#include "HMatrix.h"
#include "HFactor.h"
#include "HOutput.h"

#include <vector>
using namespace std;

class HModel {
public:
    /*
     * Tells the final solution for output
     *
     * -1 : not been set
     *  0 : completed successfully
     *  1 : infeasible
     *  2 : unbounded
     *  3 : max number of iterations
     *  4 : no solution has been found
     *  5 : max number of solutions
     *  6 : lack of storage of file spaces
     */
    int problem_status;
    double objective;

    bool solving_with_log;
    bool solving_keep_dse;
private:
    int solving_phase; // 1 = dual(1), 2 = dual(2), 4 = primal(2), other = done
    bool solving_need_invert;
    bool solving_perturbed;

public:
    HModel();
    void setup();
    void scale();
    void solve();
    void sense();
private:
    void solveInitPrimal();
    void solveInitDual();

    void solveDualp1();
    void solveDualp2();
    void solvePrimal();

    void solveDualRebuild();
    void solveDualChooseRow();
    void solveDualChooseColumn();
    void solveDualUpdate();
    void solveDualCleanup();

    void solvePrimalRebuild();
    void solvePrimalChooseColumn();
    void solvePrimalChooseRow();
    void solvePrimalUpdate();

public:
    HInput hinput;
    HOutput houtput;
private:
    // Simplex parameters
    int limit_update, limit_iteration;
    int count_update, count_iteration;
    double tol_p, tol_d;

    // LP model size, bounds, cost and scaling
    int numCol, numRow, numTotal;
    vector<double> lower0, lower1, lower_, Blower_;
    vector<double> upper0, upper1, upper_, Bupper_;
    vector<double> cost0, cost_;
    vector<double> col_scale, row_scale;

    // LP solutions
    vector<int> Nflag_, Nmove_, Bindex_;
    vector<double> dual_, value_, Bvalue_;

    // LP matrix factor and updated vectors
    HMatrix matrix;
    HFactor factor;
    HVector row_ep, row_ap, column;

    // Simplex solving buffer
    vector<int> iWork_, infeasList_;
    vector<double> dWork_, infeasArray_, infeasWeight_;

    // Simplex solving variables
    int diCount_, piCount_, rowOut_, columnOut_, columnIn_;
    double delta_, alpha_, alphaRow_, thetaDual_, thetaPrimal_;
    double row_ep_hist_dsty, column_hist_dsty, rowdse_hist_dsty;

private:
    // Random generator
    unsigned random_mw;
    unsigned random_mz;
    int intRandom() {
        random_mz = 36969 * (random_mz & 65535) + (random_mz >> 16);
        random_mw = 18000 * (random_mw & 65535) + (random_mw >> 16);
        unsigned result = (random_mz << 16) + random_mw;
        return result >> 1;
    }
    double dblRandom() {
        random_mz = 36969 * (random_mz & 65535) + (random_mz >> 16);
        random_mw = 18000 * (random_mw & 65535) + (random_mw >> 16);
        unsigned result = (random_mz << 16) + random_mw;
        return (result + 1.0) * 2.328306435454494e-10;
    }

};

#endif /* HMODEL_H_ */
