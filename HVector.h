#ifndef HVECTOR_H_
#define HVECTOR_H_

#include <vector>
using namespace std;

class HVector {
public:
    void setup(int size_);
    void clear();
    int size;
    int count;   // count of non zeros
    vector<int> index;   // index of non zeros
    vector<double> array;   // array

    void pack();
    bool packFlag;   // pack flag: do pack or not
    int packCount;   // pack count
    vector<int> packIndex;   // pack index
    vector<double> packValue;   // pack value

    vector<int> iwork;   // integer working buffer

};

#endif /* HVECTOR_H_ */
