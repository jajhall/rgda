#ifndef HINPUT_H_
#define HINPUT_H_

#include <vector>
using namespace std;

class HInput {
public:
    void loadMps(const char *filename);
    int numCol;
    int numRow;
    int AcountX;
    vector<int> Astart;
    vector<int> Aindex;
    vector<double> Avalue;
    vector<double> colLower;
    vector<double> colUpper;
    vector<double> rowLower;
    vector<double> rowUpper;
    vector<double> cost;
};

#endif /* HINPUT_H_ */
