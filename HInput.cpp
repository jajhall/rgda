#include "HInput.h"
#include "HConst.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <map>
using namespace std;

namespace {

bool HLocal_loadMPSLine(FILE *file,
                        int lmax,
                        char *line,
                        char *flag,
                        double *data) {
    int F1 = 1, F2 = 4, F3 = 14, F4 = 24, F5 = 39, F6 = 49;

    // check the buffer
    if (flag[1]) {
        flag[1] = 0;
        memcpy(&data[2], &line[F5], 8);
        data[0] = atof(&line[F6]);
        return true;
    }

    // try to read some to the line
    for (;;) {
        // Line input
        char *result = fgets(line, lmax, file);
        if (result == 0) {
            // ERROR
        }

        // Line trim   -- to delete tailing white spaces
        int lcnt = strlen(line) - 1;
        while (isspace(line[lcnt]) && lcnt >= 0)
            lcnt--; // Trim
        if (lcnt <= 0 || line[0] == '*')
            continue; // It is comment or empty..

        // Line expand -- to get data easier
        lcnt++;
        while (lcnt < F4)
            line[lcnt++] = ' '; // For row and bound row name
        if (lcnt == F4)
            line[lcnt++] = '0'; // For bound value
        line[lcnt] = '\0';

        // Done with section symbol
        if (line[0] != ' ') {
            flag[0] = line[0];
            return false;
        }

        // Read major symbol & name
        flag[0] = line[F1 + 1] == ' ' ? line[F1] : line[F1 + 1];
        memcpy(&data[1], &line[F2], 8);

        // Read 1st minor name & value to output
        memcpy(&data[2], &line[F3], 8);
        data[0] = atof(&line[F4]);

        // Keep 2nd minor name & value for future
        if (lcnt > F5)
            flag[1] = 1;
        break;
    }

    return true;
}

}

void HInput::loadMps(const char *filename) {
    // MPS file buffer
    FILE *file = fopen(filename, "r");
    const int lmax = 128;
    char line[lmax];
    char flag[2] = { 0, 0 };
    double data[3];

    // Initialisation
    numRow = 0;
    numCol = 0;
    AcountX = 0;
    double obj0 = 0;

    // Load NAME and ROWS
    HLocal_loadMPSLine(file, lmax, line, flag, data);
    HLocal_loadMPSLine(file, lmax, line, flag, data);

    vector<char> rowType;
    map<double, int> rowIndex;
    double objName = 0;
    while (HLocal_loadMPSLine(file, lmax, line, flag, data)) {
        if (flag[0] == 'N' && objName == 0) {
            objName = data[1];
        } else {
            rowType.push_back(flag[0]);
            rowIndex[data[1]] = numRow++;
        }
    }

    // Load COLUMNS
    map<double, int> colIndex;
    double lastName = 0;
    while (HLocal_loadMPSLine(file, lmax, line, flag, data)) {
        if (lastName != data[1]) { // New column
            lastName = data[1];
            colIndex[data[1]] = numCol++;
            cost.push_back(0);
            Astart.push_back(Aindex.size());
        }
        if (data[2] == objName) // Cost
            cost.back() = data[0];
        else if (data[0] != 0) {
            Aindex.push_back(rowIndex[data[2]]);
            Avalue.push_back(data[0]);
        }
    }
    Astart.push_back(Aindex.size());
    AcountX = Aindex.size();

    // Load RHS
    vector<double> RHS(numRow, 0);
    while (HLocal_loadMPSLine(file, lmax, line, flag, data)) {
        if (data[2] != objName) {
            int iRow = rowIndex[data[2]];
            RHS[iRow] = data[0];
        } else {
            obj0 = data[0]; // Objective offset
        }
    }

    // Load RANGES
    rowLower.resize(numRow);
    rowUpper.resize(numRow);
    if (flag[0] == 'R') {
        while (HLocal_loadMPSLine(file, lmax, line, flag, data)) {
            int iRow = rowIndex[data[2]];
            if (rowType[iRow] == 'L' || (rowType[iRow] == 'E' && data[0] < 0))
                rowLower[iRow] = RHS[iRow] - fabs(data[0]), rowUpper[iRow] =
                        RHS[iRow];
            else
                rowUpper[iRow] = RHS[iRow] + fabs(data[0]), rowLower[iRow] =
                        RHS[iRow];
            rowType[iRow] = 'X';
        }
    }

    // Setup bounds for row without 'RANGE'
    for (int iRow = 0; iRow < numRow; iRow++) {
        switch (rowType[iRow]) {
        case 'L':
            rowLower[iRow] = -H_INF, rowUpper[iRow] = RHS[iRow];
            break;
        case 'G':
            rowLower[iRow] = RHS[iRow], rowUpper[iRow] = +H_INF;
            break;
        case 'E':
            rowLower[iRow] = RHS[iRow], rowUpper[iRow] = RHS[iRow];
            break;
        case 'N':
            rowLower[iRow] = -H_INF, rowUpper[iRow] = +H_INF;
            break;
        case 'X':
            break;
        }
    }

    // Load BOUNDS
    colLower.assign(numCol, 0);
    colUpper.assign(numCol, H_INF);
    if (flag[0] == 'B') {
        while (HLocal_loadMPSLine(file, lmax, line, flag, data)) {
            int iCol = colIndex[data[2]];

            switch (flag[0]) {
            case 'O': /*LO*/
                colLower[iCol] = data[0];
                break;
            case 'I': /*MI*/
                colLower[iCol] = -H_INF;
                break;
            case 'L': /*PL*/
                colUpper[iCol] = H_INF;
                break;
            case 'X': /*FX*/
                colLower[iCol] = data[0], colUpper[iCol] = data[0];
                break;
            case 'R': /*FR*/
                colLower[iCol] = -H_INF, colUpper[iCol] = H_INF;
                break;
            case 'P': /*UP*/
                colUpper[iCol] = data[0];
                if (colLower[iCol] == 0 && data[0] < 0)
                    colLower[iCol] = -H_INF;
                break;
            }
        }
    }

    // Load ENDATA and close file
    fclose(file);
}
