#ifndef HOUTPUT_H_
#define HOUTPUT_H_

#include <vector>
using namespace std;

class HOutput {
public:
    // LP model size
    int numCol;
    int numRow;
    int numTotal;

    // LP results
    int problem_status;
    double objective;

    // Primal and dual solution
    vector<double> value;
    vector<double> dual;

    vector<int> Nflag; // Non basic or not
    vector<int> Nmove; // Non basic moving direction

    // SOBJ data
    vector<double> c_up_c;
    vector<double> c_up_f;
    vector<int> c_up_e;
    vector<int> c_up_l;
    vector<double> c_dn_c;
    vector<double> c_dn_f;
    vector<int> c_dn_e;
    vector<int> c_dn_l;

    // SBND data
    vector<double> b_up_b;
    vector<double> b_up_f;
    vector<int> b_up_e;
    vector<int> b_up_l;
    vector<double> b_dn_b;
    vector<double> b_dn_f;
    vector<int> b_dn_e;
    vector<int> b_dn_l;
};

#endif /* HOUTPUT_H_ */
