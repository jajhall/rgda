#include "HVector.h"

void HVector::setup(int size_) {
    size = size_;
    count = 0;
    index.resize(size);
    array.assign(size, 0);
    iwork.assign(size * 4, 0);

    packCount = 0;
    packIndex.resize(size);
    packValue.resize(size);
}

void HVector::clear() {
    if (count * 3 > size) {
        array.assign(size, 0);
    } else {
        for (int i = 0; i < count; i++)
            array[index[i]] = 0;
    }
    packFlag = false;
    count = 0;
}

void HVector::pack() {
    if (packFlag) {
        packFlag = false;
        packCount = count;
        for (int i = 0; i < count; i++) {
            const int ipack = index[i];
            packIndex[i] = ipack;
            packValue[i] = array[ipack];
        }
    }
}

