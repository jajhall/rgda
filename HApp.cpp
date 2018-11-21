#include "HModel.h"
#include "HInput.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
using namespace std;

extern HModel hmodel;

#include "hems.h"

void checkSolve(HModel *model) {
  model->setup();
  //  model->scale();
  model->solve();
}

void test_mps(char *filename) {
        const double HEMS_INF = 1e31;
	// Input all data in MPS format
	HInput mps_input;

	mps_input.loadMps(filename);
	int numCol = mps_input.numCol;
	int numRow = mps_input.numRow;
	int AcountX = mps_input.AcountX;
	double *cost = &mps_input.cost[0];
	double *rowLower = &mps_input.rowLower[0];
	double *rowUpper = &mps_input.rowUpper[0];
	double *colLower = &mps_input.colLower[0];
	double *colUpper = &mps_input.colUpper[0];
	int *Aindex = &mps_input.Aindex[0];
	int *Astart = &mps_input.Astart[0];
	double *Avalue = &mps_input.Avalue[0];
	for (int i = 0; i <= numCol; i++)
		Astart[i]++;
	for (int k = 0; k < AcountX; k++)
		Aindex[k]++;

	// Load, scale and solve
	int returnCode = 0;
	hems_lmdl(&returnCode, 0, 2, numRow, numCol, AcountX, &cost[0],
			&rowLower[0], &rowUpper[0], &colLower[0], &colUpper[0], &Aindex[0],
			&Astart[0], &Avalue[0]);
	//	hems_scal(&returnCode, 0);
	hems_sslv(&returnCode, 0, 0, 0);

	for (int i = 0; i <= numCol; i++)
		Astart[i]--;
	for (int k = 0; k < AcountX; k++)
		Aindex[k]--;

	// Now we can get the solution
	vector<double> colValue(numCol), colDual(numCol);
	vector<double> rowValue(numRow), rowDual(numRow);
	hems_gepr(&returnCode, 0, numRow + 1, numCol + 1, 0, 0, &rowValue[0],
			&colValue[0]);
	hems_gedu(&returnCode, 0, numRow + 1, numCol + 1, 0, 0, &rowDual[0],
			&colDual[0]);

	// Now do the ranging and get result
	vector<double> co_up_c(numCol), co_dn_c(numCol);
	vector<double> co_up_f(numCol), co_dn_f(numCol);
	vector<int> co_up_e(numCol), co_dn_e(numCol);
	vector<int> co_up_l(numCol), co_dn_l(numCol);

	vector<double> cb_up_b(numCol), cb_dn_b(numCol);
	vector<double> cb_up_f(numCol), cb_dn_f(numCol);
	vector<int> cb_up_e(numCol), cb_dn_e(numCol);
	vector<int> cb_up_l(numCol), cb_dn_l(numCol);

	vector<double> rb_up_b(numRow), rb_dn_b(numRow);
	vector<double> rb_up_f(numRow), rb_dn_f(numRow);
	vector<int> rb_up_e(numRow), rb_dn_e(numRow);
	vector<int> rb_up_l(numRow), rb_dn_l(numRow);

	hems_rgda(&returnCode, 0);
	hems_gerg(&returnCode, 0, numRow + 1, numCol + 1, 0, 0, &co_up_c[0],
			&co_dn_c[0], &co_up_f[0], &co_dn_f[0], &co_up_e[0], &co_dn_e[0],
			&co_up_l[0], &co_dn_l[0], &cb_up_b[0], &cb_dn_b[0], &cb_up_f[0],
			&cb_dn_f[0], &cb_up_e[0], &cb_dn_e[0], &cb_up_l[0], &cb_dn_l[0],
			&rb_up_b[0], &rb_dn_b[0], &rb_up_f[0], &rb_dn_f[0], &rb_up_e[0],
			&rb_dn_e[0], &rb_up_l[0], &rb_dn_l[0]);

	double totalSensitivityDataError = 0;
	const double toleranceRelativeTotalError = 1e-15;
	const double toleranceRelativeError = toleranceRelativeTotalError;
	const double relativeErrorDenominator = max(1.0, abs(hmodel.objective));
	bool reportSensitivityDataCheck = numCol+numRow < 100;

	vector<int> Nflag = hmodel.houtput.Nflag;
	vector<int> Nmove = hmodel.houtput.Nmove;

	// Show all rowwise data
	printf(" --- Row bounds ranging ---\n");
	printf("Row %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
			"lower", "upper", "value", "cost", "dual", "bound^", "object^",
			"verify^", "bound_", "object_", "verify_");
	for (int i = 0; i < numRow; i++) {
		double solved_up = 0;
		double solved_dn = 0;
		{
			if (rb_dn_b[i] > -HEMS_INF) {
			  HModel tmodel;
			  tmodel.solving_with_log = false;
			  tmodel.hinput = mps_input;
			  if (Nflag[i + numCol]) {
			    if (Nmove[i + numCol] == 0) {
			      tmodel.hinput.rowLower[i] = rb_dn_b[i];
			      tmodel.hinput.rowUpper[i] = rb_dn_b[i];
			    } else if (Nmove[i + numCol] == -1) {
			      tmodel.hinput.rowLower[i] = rb_dn_b[i];
			    } else {
			      tmodel.hinput.rowUpper[i] = rb_dn_b[i];
			    }
			  } else {
			    tmodel.hinput.rowUpper[i] = rb_dn_b[i];
			  }
			  checkSolve(&tmodel);
			  solved_dn = tmodel.objective;
			} else {
			  solved_dn = rb_dn_f[i];
			}
			double error = abs(solved_dn-rb_dn_f[i]);
			totalSensitivityDataError += error;
			double relativeError = error/relativeErrorDenominator;
			if (relativeError > toleranceRelativeError) {
			  printf("Row bound down ranging error: Row %2d; solved_dn = %12g; b_dn_f = %12g; relativeError = %12g\n",
				 i, solved_dn, rb_dn_f[i], relativeError);
			}
		}

		{

			if (rb_up_b[i] < HEMS_INF) {
			  HModel tmodel;
			  tmodel.solving_with_log = false;
			  tmodel.hinput = mps_input;
			  if (Nflag[i + numCol]) {
			    if (Nmove[i + numCol] == 0) {
			      tmodel.hinput.rowLower[i] = rb_up_b[i];
			      tmodel.hinput.rowUpper[i] = rb_up_b[i];
			    } else if (Nmove[i + numCol] == -1) {
			      tmodel.hinput.rowLower[i] = rb_up_b[i];
			    } else {
			      tmodel.hinput.rowUpper[i] = rb_up_b[i];
			    }
			  } else {
			    tmodel.hinput.rowLower[i] = rb_up_b[i];
			  }
			  checkSolve(&tmodel);
			  solved_up = tmodel.objective;
			} else {
			  solved_up = rb_up_f[i];
			}
			double error = abs(solved_up-rb_up_f[i]);
			totalSensitivityDataError += error;
			double relativeError = error/relativeErrorDenominator;
			if (relativeError > toleranceRelativeError) {
			  printf("Row bound   up ranging error: Row %2d; solved_up = %12g; b_up_f = %12g; relativeError = %12g\n",
				 i, solved_up, rb_up_f[i], relativeError);
			}
		}
		printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
				i, rowLower[i], rowUpper[i], rowValue[i], 0.0, rowDual[i],
				rb_up_b[i], rb_up_f[i], solved_up, rb_dn_b[i], rb_dn_f[i],
				solved_dn);
	}
	printf("\n\n");

	printf(" --- Column bounds ranging ---\n");
	printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
			"lower", "upper", "value", "cost", "dual", "bound^", "object^",
			"verify^", "bound_", "object_", "verify_");
	for (int i = 0; i < numCol; i++) {
		double solved_up = 0;
		double solved_dn = 0;
		{
    		        if (cb_dn_b[i] > -HEMS_INF) {
			  HModel tmodel;
			  tmodel.solving_with_log = false;
			  tmodel.hinput = mps_input;
			  if (Nflag[i]) {
			    if (Nmove[i] == 0) {
			      tmodel.hinput.colLower[i] = cb_dn_b[i];
			      tmodel.hinput.colUpper[i] = cb_dn_b[i];
			    } else if (Nmove[i] == 1) {
			      tmodel.hinput.colLower[i] = cb_dn_b[i];
			    } else {
			      tmodel.hinput.colUpper[i] = cb_dn_b[i];
			    }
			  } else {
			    tmodel.hinput.colUpper[i] = cb_dn_b[i];
			  }
			  checkSolve(&tmodel);
			  solved_dn = tmodel.objective;
			} else {
			  solved_dn = cb_dn_f[i];
			}
			double error = abs(solved_dn-cb_dn_f[i]);
			totalSensitivityDataError += error;
			double relativeError = error/relativeErrorDenominator;
			if (relativeError > toleranceRelativeError) {
			  printf("Col bound down ranging error: Col %2d; solved_dn = %12g; b_dn_f = %12g; relativeError = %12g\n",
				 i, solved_dn, cb_dn_f[i], relativeError);
			}
//#endif
		}

		{

 		        if (cb_up_b[i] < HEMS_INF) {
			  HModel tmodel;
			  tmodel.solving_with_log = false;
			  tmodel.hinput = mps_input;
			  if (Nflag[i]) {
			    if (Nmove[i] == 0) {
			      tmodel.hinput.colLower[i] = cb_up_b[i];
			      tmodel.hinput.colUpper[i] = cb_up_b[i];
			    } else if (Nmove[i] == 1) {
			      tmodel.hinput.colLower[i] = cb_up_b[i];
			    } else {
			      tmodel.hinput.colUpper[i] = cb_up_b[i];
			    }
			  } else {
			    tmodel.hinput.colLower[i] = cb_up_b[i];
			  }
			  checkSolve(&tmodel);
			  solved_up = tmodel.objective;
			} else {
			  solved_up = cb_up_f[i];
			}
			double error = abs(solved_up-cb_up_f[i]);
			totalSensitivityDataError += error;
			double relativeError = error/relativeErrorDenominator;
			if (relativeError > toleranceRelativeError) {
			  printf("Col bound   up ranging error: Col %2d; solved_up = %12g; b_up_f = %12g; relativeError = %12g\n",
				 i, solved_up, cb_up_f[i], relativeError);
			}
		}
		printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
				i, colLower[i], colUpper[i], colValue[i], cost[i], colDual[i],
				cb_up_b[i], cb_up_f[i], solved_up, cb_dn_b[i], cb_dn_f[i],
				solved_dn);
	}
	printf("\n\n");

	printf("--- Column cost ranging ---\n");

	printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
			"lower", "upper", "value", "cost", "dual", "cost^", "object^",
			"verify^", "cost_", "object_", "verify_");
	for (int i = 0; i < numCol; i++) {
		double solved_up = 0;
		double solved_dn = 0;
		{
		        if (fabs(co_dn_c[i]) < 1e30) {
			  HModel tmodel;
			  tmodel.solving_with_log = false;
			  tmodel.hinput = mps_input;
			  tmodel.hinput.cost[i] = co_dn_c[i];
			  checkSolve(&tmodel);
			  solved_dn = tmodel.objective;
			} else {
			  solved_dn = co_dn_f[i];
			}
			double error = abs(solved_dn-co_dn_f[i]);
			totalSensitivityDataError += error;
			double relativeError = error/relativeErrorDenominator;
			if (relativeError > toleranceRelativeError) {
			  printf("Col  cost down ranging error: Col %2d; solved_dn = %12g; c_dn_f = %12g; relativeError = %12g\n",
				 i, solved_dn, co_dn_f[i], relativeError);
			}
		}

		{

		        if (fabs(co_up_c[i]) < 1e30) {
			  HModel tmodel;
			  tmodel.solving_with_log = false;
			  tmodel.hinput = mps_input;
			  tmodel.hinput.cost[i] = co_up_c[i];
			  checkSolve(&tmodel);
			  solved_up = tmodel.objective;
			} else {
			  solved_up = co_up_f[i];
			}
			double error = abs(solved_up-co_up_f[i]);
			totalSensitivityDataError += error;
			double relativeError = error/relativeErrorDenominator;
			if (relativeError > toleranceRelativeError) {
			  printf("Col  cost   up ranging error: Col %2d; solved_up = %12g; c_up_f = %12g; relativeError = %12g\n",
				 i, solved_up, co_up_f[i], relativeError);
			}
		}
		printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
				i, colLower[i], colUpper[i], colValue[i], cost[i], colDual[i],
				co_up_c[i], co_up_f[i], solved_up, co_dn_c[i], co_dn_f[i],
				solved_dn);
	}
	printf("\n\n");
	totalSensitivityDataError /= relativeErrorDenominator;
	if (totalSensitivityDataError > toleranceRelativeTotalError) {
	  printf("WARNING totalSensitivityDataError = %12g\n", totalSensitivityDataError);
	}

	cout << endl;
}

int main(int argc, char **argv) {
	test_mps(argv[1]);
	return 0;
}

