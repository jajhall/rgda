#include "HModel.h"
#include "HConst.h"

#include <cstdio>
#include <cmath>
#include <set>
using namespace std;

HModel::HModel() {
	// Setup random generator
	random_mw = 1985;
	random_mz = 2012;

	// Setup parameters
	limit_update = 100;
	limit_iteration = 100000000;

	// Setup flags
	problem_status = -1;
	solving_with_log = true;
	solving_keep_dse = true;
	solving_perturbed = false;

	// Setup tolerances
	tol_p = 1e-7;
	tol_d = 1e-7;
}

void HModel::setup() {
	count_iteration = 0;

	// Copy model size
	numCol = hinput.numCol;
	numRow = hinput.numRow;
	numTotal = numCol + numRow;

	// Setup matrix
	int *Astart = &hinput.Astart[0];
	int *Aindex = &hinput.Aindex[0];
	double *Avalue = &hinput.Avalue[0];
	matrix.setup(numCol, numRow, Astart, Aindex, Avalue);

	// Setup bounds
	lower0.resize(numTotal);
	upper0.resize(numTotal);
	cost0.resize(numTotal);
	for (int i = 0; i < numCol; i++) {
		lower0[i] = hinput.colLower[i];
		upper0[i] = hinput.colUpper[i];
		cost0[i] = hinput.cost[i];
	}
	for (int i = 0, j = numCol; i < numRow; i++, j++) {
		lower0[j] = -hinput.rowUpper[i];
		upper0[j] = -hinput.rowLower[i];
		cost0[j] = 0;
	}
	lower_ = lower1 = lower0;
	upper_ = upper1 = upper0;
	dual_ = cost_ = cost0;

	for (int i = 0; i < numTotal; i++) {
		double lower = lower1[i], upper = upper1[i];
		if (lower == -H_INF && upper == H_INF) {
			lower1[i] = -1000, upper1[i] = 1000; // FREE
		} else if (lower == -H_INF) {
			lower1[i] = -1, upper1[i] = 0; // UPPER
		} else if (upper == H_INF) {
			lower1[i] = 0, upper1[i] = 1; // LOWER
		} else {
			lower1[i] = 0, upper1[i] = 0; // BOXED or FIXED
		}
	}

	// Setup scaling
	col_scale.assign(numCol, 1);
	row_scale.assign(numRow, 1);

	// Setup solution spaces
	Nflag_.assign(numTotal, 0);
	Nmove_.resize(numTotal);
	value_.resize(numTotal);
	for (int i = 0; i < numCol; i++)
		Nflag_[i] = 1;

	Bindex_.resize(numRow);
	Bvalue_.resize(numRow);
	Blower_.resize(numRow);
	Bupper_.resize(numRow);
	for (int iRow = 0; iRow < numRow; iRow++)
		Bindex_[iRow] = iRow + numCol;

	// Setup solving buffers
	iWork_.resize(8 * numTotal);
	dWork_.resize(8 * numTotal);

	infeasList_.resize(numRow);
	infeasArray_.resize(numRow);
	infeasWeight_.assign(numRow, 1);

	column.setup(numRow);
	row_ep.setup(numRow);
	row_ap.setup(numCol);
	column_hist_dsty = 0;
	row_ep_hist_dsty = 0;
	rowdse_hist_dsty = 0;

	// Setup factor
	factor.setup(&matrix);
	limit_update = min(100 + numRow / 200, 1000);
}

void HModel::scale() {
	// Reset all scaling to 1
	col_scale.assign(numCol, 1);
	row_scale.assign(numRow, 1);

	// Alias to the matrix
	const int *Astart = &matrix.Astart_[0];
	const int *Aindex = &matrix.Aindex_[0];
	const double *Avalue = &matrix.Avalue_[0];

	// Find out min0 / max0, skip on if in [0.2, 5]
	double min0 = H_INF, max0 = 0;
	for (int k = 0, AnX = Astart[numCol]; k < AnX; k++) {
		double value = fabs(Avalue[k]);
		min0 = min(min0, value);
		max0 = max(max0, value);
	}
	if (min0 >= 0.2 && max0 <= 5)
		return;

	// See if we want to include cost include if min-cost < 0.1
	double minc = H_INF;
	for (int i = 0; i < numCol; i++)
		if (cost_[i])
			minc = min(minc, fabs(cost_[i]));
	bool doCost = minc < 0.1;

	// Search up to 6 times
	vector<double> rowMin(numRow, H_INF);
	vector<double> rowMax(numRow, 1 / H_INF);
	for (int search_count = 0; search_count < 6; search_count++) {
		// Find column scale, prepare row data
		for (int iCol = 0; iCol < numCol; iCol++) {
			// For column scale (find)
			double colMin = H_INF;
			double colMax = 1 / H_INF;
			double myCost = fabs(cost_[iCol]);
			if (doCost && myCost != 0)
				colMin = min(colMin, myCost), colMax = max(colMax, myCost);
			for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
				double value = fabs(Avalue[k]) * row_scale[Aindex[k]];
				colMin = min(colMin, value), colMax = max(colMax, value);
			}
			col_scale[iCol] = 1 / sqrt(colMin * colMax);

			// For row scale (only collect)
			for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
				int iRow = Aindex[k];
				double value = fabs(Avalue[k]) * col_scale[iCol];
				rowMin[iRow] = min(rowMin[iRow], value);
				rowMax[iRow] = max(rowMax[iRow], value);
			}
		}

		// For row scale (find)
		for (int iRow = 0; iRow < numRow; iRow++) {
			row_scale[iRow] = 1 / sqrt(rowMin[iRow] * rowMax[iRow]);
		}
		rowMin.assign(numRow, H_INF);
		rowMax.assign(numRow, 1 / H_INF);
	}

	// Make it numerical better
	const double ln2 = log(2.0);
	for (int iCol = 0; iCol < numCol; iCol++)
		col_scale[iCol] = pow(2.0, floor(log(col_scale[iCol]) / ln2 + 0.5));
	for (int iRow = 0; iRow < numRow; iRow++)
		row_scale[iRow] = pow(2.0, floor(log(row_scale[iRow]) / ln2 + 0.5));

	// Apply scaling to matrix and bounds
	matrix.scale(&col_scale[0], &row_scale[0]);
	for (int iCol = 0; iCol < numCol; iCol++) {
		lower0[iCol] /= lower0[iCol] == -H_INF ? 1 : col_scale[iCol];
		upper0[iCol] /= upper0[iCol] == +H_INF ? 1 : col_scale[iCol];
		cost0[iCol] *= col_scale[iCol];
	}
	for (int iRow = 0, iTot = numCol; iRow < numRow; iRow++, iTot++) {
		lower0[iTot] *= lower0[iTot] == -H_INF ? 1 : row_scale[iRow];
		upper0[iTot] *= upper0[iTot] == +H_INF ? 1 : row_scale[iRow];
	}

	// Update working bounds
	lower_ = lower0;
	upper_ = upper0;
	cost_ = cost0;
}

void HModel::solve() {
	solveInitPrimal();
	solveInitDual();

	solving_phase = diCount_ > 0 ? 1 : 2;
	while (solving_phase) {
		switch (solving_phase) {
		case 1:
			solveDualp1();
			break;
		case 2:
			solveDualp2();
			break;
		case 4:
			solvePrimal();
			break;
		default:
			solving_phase = 0;
			break;
		}
	}

	// Fill in output buffer
	houtput.numCol = numCol;
	houtput.numRow = numRow;
	houtput.numTotal = numTotal;

	// Take primal solution
	houtput.value = value_;
	for (int iRow = 0; iRow < numRow; iRow++)
		houtput.value[Bindex_[iRow]] = Bvalue_[iRow];

	// Take dual solution
	houtput.dual = dual_;
	for (int iRow = 0; iRow < numRow; iRow++)
		houtput.dual[Bindex_[iRow]] = 0;

	// Take non basic flag and move
	houtput.Nflag = Nflag_;
	houtput.Nmove = Nmove_;

	// Scale back
	for (int iCol = 0; iCol < numCol; iCol++) {
		houtput.value[iCol] *= col_scale[iCol];
		houtput.dual[iCol] /= col_scale[iCol];

	}
	for (int iRow = 0, iTot = numCol; iRow < numRow; iRow++, iTot++) {
		houtput.value[iTot] /= row_scale[iRow];
		houtput.dual[iTot] *= row_scale[iRow];
	}

}

void HModel::solveInitPrimal() {
	// Init primal solution
	for (int i = 0; i < numTotal; i++) {
		if (Nflag_[i]) {
			const double lower = lower_[i];
			const double upper = upper_[i];
			const double cost = dual_[i];
			if (lower == -H_INF && upper == H_INF) {
				value_[i] = 0, Nmove_[i] = 0; // Free
			} else if (lower == -H_INF) {
				value_[i] = upper, Nmove_[i] = -1; // Upper
			} else if (upper == H_INF) {
				value_[i] = lower, Nmove_[i] = 1; // Lower
			} else if (lower == upper) {
				value_[i] = lower, Nmove_[i] = 0; // Fixed
			} else if (cost >= 0) {
				value_[i] = lower, Nmove_[i] = 1; // Boxed at Lower
			} else {
				value_[i] = upper, Nmove_[i] = -1; // Boxed at Upper
			}
		} else {
			Nmove_[i] = 0;
		}
	}
}

void HModel::solveInitDual() {
	// Init dual perturbation
	set<double> cost_set;
	for (int i = 0; i < numCol; i++)
		cost_set.insert(fabs(cost_[i]));
	bool doPerturb = (int) (cost_set.size() * 4) < numCol;
	if (doPerturb) {
		solving_perturbed = true;
		double bigc = max(0.01 * (*cost_set.rbegin()), 1.0);
		double base = tol_d * 5;
		for (int i = 0; i < numCol; i++) {
			double lower = lower_[i];
			double upper = upper_[i];
			double xpert = (fabs(cost_[i]) + bigc) * base * (1 + dblRandom());
			if (lower == -H_INF && upper == H_INF) {
				// Free - no perturb
			} else if (upper == H_INF) { // Lower
				cost_[i] += xpert;
			} else if (lower == -H_INF) { // Upper
				cost_[i] += -xpert;
			} else if (lower != upper) { // Boxed
				cost_[i] += (cost_[i] >= 0) ? xpert : -xpert;
			} else {
				// Fixed - no perturb
			}
		}
		dual_ = cost_;
	}

	// Provide an initial dual infeasible count
	diCount_ = 0;
	for (int i = 0; i < numTotal; i++) {
		if (Nflag_[i]) {
			// Free
			if (lower_[i] == -H_INF && upper_[i] == H_INF) {
				diCount_ += (fabs(dual_[i]) > tol_d);
			}
			// Lower or upper
			if (lower_[i] == -H_INF || upper_[i] == H_INF) {
				diCount_ += (Nmove_[i] * dual_[i] < -tol_d);
			}
		}
	}
}

void HModel::solveDualp1() {
	if (solving_with_log)
		printf("dual-phase-1-start\n");

	// Switch to dual phase 1 bounds
	lower_ = lower1;
	upper_ = upper1;
	solveInitPrimal();

	// Main solving structure
	for (;;) {
		solveDualRebuild();
		if (solving_with_log)
			printf("%10d  %20.10e\n", count_iteration, objective);

		for (;;) {
			solveDualChooseRow();
			if (rowOut_ == -1)
				break;

			solveDualChooseColumn();
			if (columnIn_ == -1)
				break;

			solveDualUpdate();
			if (solving_need_invert)
				break;
		}

		// Fresh factorization required
		if (count_update == 0)
			break;
	}

	if (rowOut_ == -1) {
		if (solving_with_log)
			printf("dual-phase-1-optimal\n");
		// Go to phase 2
		if (objective == 0) {
			solving_phase = 2;
		} else {
			// We still have dual infeasible
			if (solving_perturbed) {
				// Clean up perturbation and go on
				solveDualCleanup();
			} else {
				// Report dual infeasible
				if (solving_with_log)
					printf("dual-infeasible\n");
				solving_phase = -1;
				problem_status = 2;
			}
		}
	} else if (columnIn_ == -1) {
		if (solving_with_log)
			printf("dual-phase-1-unbounded\n");
		// We got dual phase 1 unbounded - strange
		if (solving_perturbed) {
			// Clean up perturbation and go on
			solveDualCleanup();
		} else {
			// Report strange issues
			if (solving_with_log)
				printf("dual-phase-1-not-solved\n");
			solving_phase = -1;
			problem_status = 4;
		}
	}

	if (solving_phase == 2) {
		lower_ = lower0;
		upper_ = upper0;
		solveInitPrimal();
	}

}

void HModel::solveDualp2() {
	if (solving_with_log)
		printf("dual-phase-2-start\n");

	// Main solving structure
	for (;;) {
		solveDualRebuild();
		if (diCount_ > 0)
			break;
		if (solving_with_log)
			printf("%10d  %20.10e\n", count_iteration, objective);

		for (;;) {
			solveDualChooseRow();
			if (rowOut_ == -1)
				break;

			solveDualChooseColumn();
			if (columnIn_ == -1)
				break;

			solveDualUpdate();
			if (solving_need_invert)
				break;
		}

		// Fresh factorization required
		if (count_update == 0)
			break;
	}

	if (diCount_ > 0) {
		// We got free variables - need dual phase 1
		if (solving_with_log)
			printf("dual-phase-2-found-free\n");
		solving_phase = 1;
	} else if (rowOut_ == -1) {
		// Tentative optimal
		if (solving_with_log)
			printf("dual-phase-2-optimal\n");

		// See if we need primal
		solveDualCleanup();
		if (diCount_ > 0) {
			solving_phase = 4;
		} else {
			solving_phase = 0;
			if (solving_with_log)
				printf("problem-optimal\n");
			problem_status = 0;
		}
	} else if (columnIn_ == -1) {
		if (solving_with_log)
			printf("dual-phase-2-unbounded\n");
		// We got dual phase 1 unbounded - strange
		if (solving_perturbed) {
			// Clean up perturbation and go on
			solveDualCleanup();
		} else {
			// Report strange issues
			if (solving_with_log)
				printf("problem-infeasible\n");
			solving_phase = -1;
			problem_status = 1;
		}
	}
}

void HModel::solvePrimal() {
	if (solving_with_log)
		printf("primal-start\n");

	for (;;) {
		solvePrimalRebuild();
		if (solving_with_log)
			printf("%10d  %20.10e\n", count_iteration, objective);

		for (;;) {
			solvePrimalChooseColumn();
			if (columnIn_ == -1)
				break;

			solvePrimalChooseRow();
			if (rowOut_ == -1) {
				break;
			}
			solvePrimalUpdate();
			if (solving_need_invert)
				break;
		}

		// Fresh factorization required
		if (count_update == 0)
			break;
	}
	if (columnIn_ == -1) {
		if (solving_with_log)
			printf("primal-optimal\n");

		if (solving_with_log)
			printf("problem-optimal\n");
		problem_status = 0;
	} else {
		// Need to consider shift
		if (solving_with_log)
			printf("problem-unbounded\n");
		problem_status = 2;
	}

	solving_phase = -1;
}

void HModel::solveDualRebuild() {
	solving_need_invert = false;

	// Rebuild factor
	if (count_iteration == 0 || count_update > 0) {
		count_update = 0;
		if (solving_keep_dse) {
			for (int i = 0; i < numRow; i++)
				dWork_[Bindex_[i]] = infeasWeight_[i];
		}

		factor.build(&Bindex_[0]);

		if (solving_keep_dse) {
			for (int i = 0; i < numRow; i++)
				infeasWeight_[i] = dWork_[Bindex_[i]];
		}
	}

	// Recompute dual solution
	row_ep.clear();
	row_ap.clear();
	for (int iRow = 0; iRow < numRow; iRow++) {
		row_ep.index[iRow] = iRow;
		row_ep.array[iRow] = cost_[Bindex_[iRow]];
	}
	row_ep.count = numRow;
	factor.btran(row_ep, 1);
	matrix.price_by_col(row_ap, row_ep);
	for (int i = 0; i < numCol; i++)
		dual_[i] = cost_[i] - row_ap.array[i];
	for (int i = numCol; i < numTotal; i++)
		dual_[i] = cost_[i] - row_ep.array[i - numCol];

	// Find out any infeasible (flip or shift) - only free are recorded
	diCount_ = 0;
	for (int i = 0; i < numTotal; i++) {
		if (Nflag_[i]) {
			if (lower_[i] == -H_INF && upper_[i] == H_INF) {
				// FREE variable
				diCount_ += (fabs(dual_[i]) > tol_d);
			} else if (Nmove_[i] * dual_[i] < -tol_d) {
				if (lower_[i] != -H_INF && upper_[i] != H_INF) {
					// Boxed variable = flip
					Nmove_[i] *= -1;
					value_[i] = dual_[i] > 0 ? lower_[i] : upper_[i];
				} else {
					// Other variable = shift
					solving_perturbed = true;
					if (Nmove_[i] == 1) {
						double dual = (1 + dblRandom()) * tol_d;
						double shift = dual - dual_[i];
						dual_[i] = dual;
						cost_[i] = cost_[i] + shift;
					} else {
						double dual = -(1 + dblRandom()) * tol_d;
						double shift = dual - dual_[i];
						dual_[i] = dual;
						cost_[i] = cost_[i] + shift;
					}
				}
			}
		}
	}

	// Compute basic primal solutions
	column.clear();
	for (int i = 0; i < numTotal; i++)
		if (Nflag_[i] && value_[i] != 0)
			matrix.collect_aj(column, i, value_[i]);
	factor.ftran(column, 1);
	for (int i = 0; i < numRow; i++) {
		const int iCol = Bindex_[i];
		Blower_[i] = lower_[iCol];
		Bupper_[i] = upper_[iCol];
		Bvalue_[i] = -column.array[i];
	}

	// Collect primal infeasible as a list
	piCount_ = 0;
	for (int i = 0; i < numRow; i++) {
		const double my_lower = Blower_[i];
		const double my_upper = Bupper_[i];
		const double my_value = Bvalue_[i];
		double infeas = 0;
		if (my_lower - my_value > tol_p) {
			infeas = my_lower - my_value;
			infeasList_[piCount_++] = i;
		} else if (my_value - my_upper > tol_p) {
			infeas = my_value - my_upper;
			infeasList_[piCount_++] = i;
		}
		infeasArray_[i] = infeas * infeas;
	}

	// Compute the objective value
	objective = 0;
	for (int i = 0; i < numTotal; i++)
		if (Nflag_[i])
			objective += value_[i] * dual_[i];

}

void HModel::solveDualChooseRow() {
	// Initialize
	rowOut_ = -1;
	columnOut_ = -1;
	if (piCount_ == 0)
		return;

	// Choose a row
	int randomStart = intRandom() % piCount_;
	double bestMerit = 0;
	for (int section = 0; section < 2; section++) {
		const int start = (section == 0) ? randomStart : 0;
		const int end = (section == 0) ? piCount_ : randomStart;
		for (int i = start; i < end; i++) {
			int iRow = infeasList_[i];
			if (infeasArray_[iRow] > H_TZ) {
				const double myInfeas = infeasArray_[iRow];
				const double myWeight = infeasWeight_[iRow];
				if (bestMerit * myWeight < myInfeas) {
					bestMerit = myInfeas / myWeight;
					rowOut_ = iRow;
				}
			}
		}
	}

	// Determine delta
	if (rowOut_ != -1) {
		columnOut_ = Bindex_[rowOut_];
		if (Bvalue_[rowOut_] < Blower_[rowOut_])
			delta_ = Bvalue_[rowOut_] - Blower_[rowOut_];
		else
			delta_ = Bvalue_[rowOut_] - Bupper_[rowOut_];
	}
}

void HModel::solveDualChooseColumn() {
	// Initialize
	columnIn_ = -1;

	// Compute pivot row
	row_ep.clear();
	row_ap.clear();
	row_ep.count = 1;
	row_ep.index[0] = rowOut_;
	row_ep.array[rowOut_] = 1;
	row_ep.packFlag = true;
	factor.btran(row_ep, row_ep_hist_dsty);
	matrix.price_by_row(row_ap, row_ep);
	row_ep_hist_dsty = 0.95 * row_ep_hist_dsty + 0.05 * row_ep.count / numRow;

	// Alias to updated row
	const int ap_count = row_ap.count;
	const int ep_count = row_ep.count;
	const int *ap_index = &row_ap.index[0];
	const int *ep_index = &row_ep.index[0];
	const double *ap_array = &row_ap.array[0];
	const double *ep_array = &row_ep.array[0];

	// Alias to working buffer
	int ch1n = 0;
	int ch2n = 0;
	double ch1t = H_INF;
	double ch2t = H_INF;
	int *ch1i = &iWork_[0];
	int *ch2i = &iWork_[numCol];
	int *ch3i = &iWork_[numCol * 2];
	double *ch1x = &dWork_[0];
	double *ch2x = &dWork_[numCol];
	const double *lower = &lower_[0];
	const double *upper = &upper_[0];
	const double *dual = &dual_[0];
	const int *Nmove = &Nmove_[0];
	const double totalDelta = fabs(delta_);
	const int deltaSign = delta_ < 0 ? -1 : 1;
	double totalChange = 0;

	// Choose column pass 1
	double tol_a = count_update < 10 ? 1e-9 : count_update < 20 ? 1e-8 : 1e-7;
	for (int i = 0; i < ap_count; i++) {
		int iCol = ap_index[i];
		double my_value = ap_array[iCol] * Nmove[iCol] * deltaSign;
		if (my_value > tol_a) {
			ch1i[ch1n] = iCol;
			ch1x[ch1n++] = my_value;
			double cost = dual[iCol] * Nmove[iCol];
			double compare = my_value * ch1t - cost;
			if (compare > tol_d)
				ch1t = (cost + tol_d) / my_value;
		}
	}
	for (int i = 0; i < ep_count; i++) {
		int iGet = ep_index[i];
		int iCol = iGet + numCol;
		double my_value = ep_array[iGet] * Nmove[iCol] * deltaSign;
		if (my_value > tol_a) {
			ch1i[ch1n] = iCol;
			ch1x[ch1n++] = my_value;
			double cost = dual[iCol] * Nmove[iCol];
			double compare = my_value * ch1t - cost;
			if (compare > tol_d)
				ch1t = (cost + tol_d) / my_value;
		}
	}

	if (ch1n == 0)
		return;

	// Choose column pass 2
	ch1t = max(10 * ch1t, 1e-7);
	while (ch1t < 1e18) {
		ch2n = 0;
		ch2t = H_INF;
		double *ch2x_p = ch2x;
		totalChange = 0;

		for (int i = 0; i < ch1n; i++) {
			int icol = ch1i[i];
			double alpha = ch1x[i];
			double cost = dual[icol] * Nmove[icol];
			double selectCompare = alpha * ch1t - cost;
			if (selectCompare > tol_d) { // Select
				ch2i[ch2n++] = icol;
				*ch2x_p++ = (cost + tol_d) / alpha; // Relax
				*ch2x_p++ = cost / alpha; // Tight
				*ch2x_p++ = alpha;
				*ch2x_p++ = upper[icol] - lower[icol];
				totalChange += alpha * (upper[icol] - lower[icol]);
			} else { // Update upper_theta
				double compare = alpha * ch2t - cost;
				if (compare > tol_d)
					ch2t = (cost + tol_d) / alpha;
			}
		}

		ch1t = max(2 * ch2t, 2 * ch1t);
		if (totalChange > totalDelta || ch2n == ch1n)
			break;
	}

	// Choose column pass 3
	int ch2left = ch2n;
	int *ch2mark = &ch3i[0];
	for (int i = 0; i < ch2n; i++)
		ch2mark[i] = 1;
	int ch3np = 0;
	int *ch3ip = &ch3i[ch2n];

	totalChange = 0;
	while (ch2left > 0) {
		// Get upperTheta of the available part
		double upperTheta = 1e200;
		for (int ich2 = 0; ich2 < ch2n; ich2++) {
			const double relax = ch2x[ich2 * 4];
			if (ch2mark[ich2] && relax < upperTheta)
				upperTheta = relax;
		}

		int bestip = -1;
		double bestx = 0;
		for (int ich2 = 0; ich2 < ch2n; ich2++) {
			const int iTight = ich2 * 4 + 1;
			const double tight = ch2x[iTight];
			if (ch2mark[ich2] && tight <= upperTheta) {
				ch2mark[ich2] = 0;
				ch2left--;
				const double alpha = ch2x[ich2 * 4 + 2];
				const double range = ch2x[ich2 * 4 + 3];
				totalChange += alpha * range;
				if (alpha > bestx) {
					bestx = alpha;
					bestip = ich2;
				}
			}
		}

		ch3ip[ch3np++] = bestip;
		if (totalChange >= totalDelta)
			break;
	}

	// Find max and then the tolerance
	double finalMax = 0;
	for (int iptr = 0; iptr < ch3np; iptr++) {
		int ich2 = ch3ip[iptr];
		const int iAlpha = ich2 * 4 + 2;
		if (finalMax < ch2x[iAlpha])
			finalMax = ch2x[iAlpha];
	}
	if (finalMax > 10.0)
		finalMax = 10.0;
	const double finalTol = 0.1 * finalMax;

	// Find the final choice backwardly
	for (int iptr = ch3np - 1; iptr >= 0; iptr--) {
		int ich2 = ch3ip[iptr];
		const int iAlpha = ich2 * 4 + 2;
		if (ch2x[iAlpha] > finalTol) {
			columnIn_ = ch2i[ich2];
			break;
		}
	}

	// Choose column simplex variables
	alphaRow_ = matrix.compute_dot(row_ep, columnIn_);
	if (dual_[columnIn_] * Nmove_[columnIn_] > 0) {
		thetaDual_ = dual_[columnIn_] / alphaRow_;
	} else {
		thetaDual_ = 0;
	}
}

void HModel::solveDualUpdate() {
	// Update dual
	diCount_ = 0;
	if (thetaDual_ == 0) {
		// Shift the dual to be real "EXPAND"
		solving_perturbed = true;
		cost_[columnIn_] -= dual_[columnIn_];
	} else {
		const int ap_count = row_ap.count;
		const int ep_count = row_ep.count;
		const int *ap_index = &row_ap.index[0];
		const int *ep_index = &row_ep.index[0];
		const double *ap_array = &row_ap.array[0];
		const double *ep_array = &row_ep.array[0];

		const int *Nmove = &Nmove_[0];
		double *dual = &dual_[0];
		int *diList = &iWork_[0];

		for (int i = 0; i < ap_count; i++) {
			int iCol = ap_index[i];
			dual[iCol] -= thetaDual_ * ap_array[iCol];
			if (dual[iCol] * Nmove[iCol] <= -tol_d)
				diList[diCount_++] = iCol;
		}

		Nmove = &Nmove_[numCol];
		dual = &dual_[numCol];
		for (int i = 0; i < ep_count; i++) {
			int iGet = ep_index[i];
			dual[iGet] -= thetaDual_ * ep_array[iGet];
			if (dual[iGet] * Nmove[iGet] <= -tol_d)
				diList[diCount_++] = iGet + numCol;
		}
	}
	dual_[columnIn_] = 0;
	dual_[columnOut_] = -thetaDual_;

	// Update dual - flip bounds
	column.clear();
	for (int i = 0; i < diCount_; i++) {
		int iCol = iWork_[i];
		double move = Nmove_[iCol];
		if (move * upper_[iCol] == H_INF || move * lower_[iCol] == H_INF) {
			// Can't flip to infinity, will correct by rebuild and shift
			solving_need_invert = true;
		} else {
			// Flip and collect for FTRAN-BFRT
			Nmove_[iCol] *= -1;
			value_[iCol] = Nmove_[iCol] == 1 ? lower_[iCol] : upper_[iCol];
			move *= upper_[iCol] - lower_[iCol];
			matrix.collect_aj(column, iCol, move);
			objective += move * dual_[iCol];
		}
	}

	// Update dual - BFRT to basic variables
	if (column.count > 0) {
		factor.ftran(column, column_hist_dsty);
		for (int i = 0; i < column.count; i++) {
			int iRow = column.index[i];
			Bvalue_[iRow] -= column.array[iRow];
			const double value = Bvalue_[iRow];
			const double lower = Blower_[iRow];
			const double upper = Bupper_[iRow];

			double infeas = 0;
			if (value < lower - tol_p)
				infeas = value - lower;
			if (value > upper + tol_p)
				infeas = value - upper;

			const double infeas0 = infeasArray_[iRow];
			const double infeas1 = infeas * infeas;
			if (infeas0) {
				infeasArray_[iRow] = infeas1 > 0 ? infeas1 : H_TZ;
			} else if (infeas1) {
				infeasArray_[iRow] = infeas1;
				infeasList_[piCount_++] = iRow;
			}
		}
	}

	// Update primal - basic variables
	column.clear();
	column.packFlag = true;
	matrix.collect_aj(column, columnIn_, 1);
	factor.ftran(column, column_hist_dsty);
	column_hist_dsty = 0.95 * column_hist_dsty + 0.05 * column.count / numRow;
	alpha_ = column.array[rowOut_];

	double valueOut = Bvalue_[rowOut_];
	double lowerOut = Blower_[rowOut_];
	double upperOut = Bupper_[rowOut_];
	double deltaX = delta_ < 0 ? (valueOut - lowerOut) : (valueOut - upperOut);
	thetaPrimal_ = deltaX / alpha_;
	for (int i = 0; i < column.count; i++) {
		int iRow = column.index[i];
		Bvalue_[iRow] -= thetaPrimal_ * column.array[iRow];
		const double value = Bvalue_[iRow];
		const double lower = Blower_[iRow];
		const double upper = Bupper_[iRow];

		double infeas = 0;
		if (value < lower - tol_p)
			infeas = value - lower;
		if (value > upper + tol_p)
			infeas = value - upper;

		const double infeas0 = infeasArray_[iRow];
		const double infeas1 = infeas * infeas;
		if (infeas0) {
			infeasArray_[iRow] = infeas1 > 0 ? infeas1 : H_TZ;
		} else if (infeas1) {
			infeasArray_[iRow] = infeas1;
			infeasList_[piCount_++] = iRow;
		}
	}

	// Update DSE weight - compute pivot's weight
	double pivotWeight = 0;
	for (int i = 0; i < row_ep.count; i++) {
		double value = row_ep.array[row_ep.index[i]];
		pivotWeight += value * value;
	}

	// Update DSE weight - update for others
	factor.ftran(row_ep, rowdse_hist_dsty);
	rowdse_hist_dsty = 0.95 * rowdse_hist_dsty + 0.05 * row_ep.count / numRow;
	const double weight = pivotWeight / (alpha_ * alpha_);
	const double Kai = -2 / alpha_;
	for (int i = 0; i < column.count; i++) {
		const int iRow = column.index[i];
		const double val = column.array[iRow];
		infeasWeight_[iRow] += val * (weight * val + Kai * row_ep.array[iRow]);
		if (infeasWeight_[iRow] < 1e-4)
			infeasWeight_[iRow] = 1e-4;
	}
	infeasWeight_[rowOut_] = pivotWeight / (alpha_ * alpha_);

	// Update pivot - column out
	if (lower_[columnOut_] == upper_[columnOut_]) {
		value_[columnOut_] = lower_[columnOut_];
		Nmove_[columnOut_] = 0;
	} else if (delta_ < 0) {
		value_[columnOut_] = lower_[columnOut_];
		Nmove_[columnOut_] = 1;
	} else {
		value_[columnOut_] = upper_[columnOut_];
		Nmove_[columnOut_] = -1;
	}
	Nflag_[columnOut_] = 1;

	// Update pivot - column in
	Nmove_[columnIn_] = 0;
	Nflag_[columnIn_] = 0;
	Bindex_[rowOut_] = columnIn_;
	Blower_[rowOut_] = lower_[columnIn_];
	Bupper_[rowOut_] = upper_[columnIn_];
	Bvalue_[rowOut_] = value_[columnIn_] + thetaPrimal_;
	double pivotInfeas = 0;
	if (Bvalue_[rowOut_] < Blower_[rowOut_] - tol_p)
		pivotInfeas = Bvalue_[rowOut_] - Blower_[rowOut_];
	if (Bvalue_[rowOut_] > Bupper_[rowOut_] + tol_p)
		pivotInfeas = Bvalue_[rowOut_] - Bupper_[rowOut_];
	infeasArray_[rowOut_] = max(pivotInfeas * pivotInfeas, H_TZ);

	// Update pivot - basis matrix and its factor
	factor.update(column, row_ep, rowOut_, solving_need_invert);
	matrix.update(columnIn_, columnOut_);
	if (++count_update >= limit_update)
		solving_need_invert = true;

	// Update pivot - check different alpha
	double alphaC = fabs(alpha_);
	double alphaR = fabs(alphaRow_);
	double alphaDiff = fabs(alphaC - alphaR);
	if (alphaDiff / max(alphaC, alphaR) > 1e-7)
		solving_need_invert = true;

	// Update pivot - objective
	objective += thetaDual_ * delta_;
	count_iteration++;
}

void HModel::solveDualCleanup() {
	// Remove perturbation and recompute the dual solution
	if (solving_with_log)
		printf("dual-cleanup-shift\n");
	cost_ = cost0;
	solving_perturbed = false;

	// Compute dual
	row_ep.clear();
	row_ap.clear();
	for (int iRow = 0; iRow < numRow; iRow++) {
		row_ep.index[iRow] = iRow;
		row_ep.array[iRow] = cost_[Bindex_[iRow]];
	}
	row_ep.count = numRow;
	factor.btran(row_ep, 1);
	matrix.price_by_col(row_ap, row_ep);
	for (int i = 0; i < numCol; i++)
		dual_[i] = cost_[i] - row_ap.array[i];
	for (int i = numCol; i < numTotal; i++)
		dual_[i] = cost_[i] - row_ep.array[i - numCol];

	// Count dual infeasible - no flip
	diCount_ = 0;
	for (int i = 0; i < numTotal; i++) {
		if (Nflag_[i] && fabs(dual_[i]) > tol_d) {
			diCount_ += (lower_[i] == -H_INF && upper_[i] == H_INF); // FREE
			diCount_ += (Nmove_[i] * dual_[i] < -tol_d); // Infeasible
		}
	}

	// Compute objective
	objective = 0;
	for (int i = 0; i < numTotal; i++)
		if (Nflag_[i])
			objective += value_[i] * dual_[i];

	if (solving_with_log)
		printf("%10d  %20.10e\n", count_iteration, objective);
}

void HModel::solvePrimalRebuild() {
	// Rebuild factor - only if we got updates
	solving_need_invert = false;
	if (count_iteration == 0 || count_update > 0) {
		factor.build(&Bindex_[0]);
		count_update = 0;
	}

	// Compute dual
	row_ep.clear();
	row_ap.clear();
	for (int iRow = 0; iRow < numRow; iRow++) {
		row_ep.index[iRow] = iRow;
		row_ep.array[iRow] = cost_[Bindex_[iRow]];
	}
	row_ep.count = numRow;
	factor.btran(row_ep, 1);
	matrix.price_by_col(row_ap, row_ep);
	for (int i = 0; i < numCol; i++)
		dual_[i] = cost_[i] - row_ap.array[i];
	for (int i = numCol; i < numTotal; i++)
		dual_[i] = cost_[i] - row_ep.array[i - numCol];

	// Compute basic primal
	column.clear();
	for (int i = 0; i < numTotal; i++)
		if (Nflag_[i] && value_[i] != 0)
			matrix.collect_aj(column, i, value_[i]);
	factor.ftran(column, 1);
	for (int i = 0; i < numRow; i++) {
		const int iCol = Bindex_[i];
		Blower_[i] = lower_[iCol];
		Bupper_[i] = upper_[iCol];
		Bvalue_[i] = -column.array[i];
		// TODO Check infeasible and shift

	}

	// Compute the objective value
	objective = 0;
	for (int i = 0; i < numTotal; i++)
		if (Nflag_[i])
			objective += value_[i] * dual_[i];

}

void HModel::solvePrimalChooseColumn() {
	columnIn_ = -1;
	double bestInfeas = 0;
	for (int iCol = 0; iCol < numTotal; iCol++) {
		if (Nflag_[iCol] && fabs(dual_[iCol]) > tol_d) {
			// Always take free
			// TODO: if we found free,
			// Then deal with it in dual phase 1
			if (lower_[iCol] == -H_INF && upper_[iCol] == H_INF) {
				columnIn_ = iCol;
				break;
			}
			// Then look at dual infeasible
			if (Nmove_[iCol] * dual_[iCol] < -tol_d) {
				if (bestInfeas < fabs(dual_[iCol])) {
					bestInfeas = fabs(dual_[iCol]);
					columnIn_ = iCol;
				}
			}
		}
	}

}

void HModel::solvePrimalChooseRow() {
	// Compute pivot column
	column.clear();
	column.packFlag = true;
	matrix.collect_aj(column, columnIn_, 1);
	factor.ftran(column, column_hist_dsty);
	column_hist_dsty = 0.95 * column_hist_dsty + 0.05 * column.count / numRow;

	// Initialize
	rowOut_ = -1;

	// Choose column pass 1
	double alphaTol = count_update < 10 ? 1e-9 :
						count_update < 20 ? 1e-8 : 1e-7;
	int moveIn = Nmove_[columnIn_];
	if (moveIn == 0) {
		// If there's still free in the N
		// We would report not-solved
		// Need to handle free
	}
	double relaxTheta = 1e100;
	for (int i = 0; i < column.count; i++) {
		int index = column.index[i];
		double alpha = column.array[index] * moveIn;
		if (alpha > alphaTol) {
			double relaxSpace = Bvalue_[index] - Blower_[index] + tol_p;
			if (relaxSpace < relaxTheta * alpha)
				relaxTheta = relaxSpace / alpha;
		} else if (alpha < -alphaTol) {
			double relaxSpace = Bvalue_[index] - Bupper_[index] - tol_p;
			if (relaxSpace > relaxTheta * alpha)
				relaxTheta = relaxSpace / alpha;
		}
	}

	// Choose column pass 2
	double bestAlpha = 0;
	for (int i = 0; i < column.count; i++) {
		int index = column.index[i];
		double alpha = column.array[index] * moveIn;
		if (alpha > alphaTol) {
			double tightSpace = Bvalue_[index] - Blower_[index];
			if (tightSpace < relaxTheta * alpha) {
				if (bestAlpha < alpha) {
					bestAlpha = alpha;
					rowOut_ = index;
				}
			}
		} else if (alpha < -alphaTol) {
			double tightSpace = Bvalue_[index] - Bupper_[index];
			if (tightSpace > relaxTheta * alpha) {
				if (bestAlpha < -alpha) {
					bestAlpha = -alpha;
					rowOut_ = index;
				}
			}
		}
	}

	// Other info
	if (rowOut_ != -1) {
		columnOut_ = Bindex_[rowOut_];
		alpha_ = column.array[rowOut_];
		if (alpha_ * moveIn > 0) {
			// Lower bound
			thetaPrimal_ = (Bvalue_[rowOut_] - Blower_[rowOut_]) / alpha_;
		} else {
			// Upper bound
			thetaPrimal_ = (Bvalue_[rowOut_] - Bupper_[rowOut_]) / alpha_;
		}
	}
}

void HModel::solvePrimalUpdate() {
	int moveIn = Nmove_[columnIn_];

	// 1. Make sure it is inside bounds or just flip bound
	double lowerIn = lower_[columnIn_];
	double upperIn = upper_[columnIn_];
	double valueIn = value_[columnIn_] + thetaPrimal_;
	bool flipped = false;
	if (Nmove_[columnIn_] == 1) {
		if (valueIn > upperIn + tol_p) {
			// Flip to upper
			value_[columnIn_] = upperIn;
			thetaPrimal_ = upperIn - lowerIn;
			flipped = true;
			Nmove_[columnIn_] = -1;
		}
	} else if (Nmove_[columnIn_] == -1) {
		if (valueIn < lowerIn - tol_p) {
			// Flip to lower
			value_[columnIn_] = lowerIn;
			thetaPrimal_ = lowerIn - upperIn;
			flipped = true;
			Nmove_[columnIn_] = 1;
		}
	}

	for (int i = 0; i < column.count; i++) {
		int index = column.index[i];
		Bvalue_[index] -= thetaPrimal_ * column.array[index];
	}

	// If flipped, then no need touch the pivots
	if (flipped) {
		count_iteration++;
		objective = 0;
		for (int i = 0; i < numTotal; i++)
			if (Nflag_[i])
				objective += value_[i] * dual_[i];
		return;
	}

	// Pivot in
	Bindex_[rowOut_] = columnIn_;
	Blower_[rowOut_] = lowerIn;
	Bupper_[rowOut_] = upperIn;
	Bvalue_[rowOut_] = valueIn;
	Nflag_[columnIn_] = 0;
	Nmove_[columnIn_] = 0;

	// Check for any possible infeasible
	for (int iRow = 0; iRow < numRow; iRow++) {
		if (Bvalue_[iRow] < Blower_[iRow] - tol_p) {
			solving_need_invert = true;
		} else if (Bvalue_[iRow] > Bupper_[iRow] + tol_p) {
			solving_need_invert = true;
		}
	}

	// Pivot out
	if (lower_[columnOut_] == upper_[columnOut_]) {
		value_[columnOut_] = lower_[columnOut_];
		Nmove_[columnOut_] = 0;
	} else if (alpha_ * moveIn > 0) {
		value_[columnOut_] = lower_[columnOut_];
		Nmove_[columnOut_] = 1;
	} else {
		value_[columnOut_] = upper_[columnOut_];
		Nmove_[columnOut_] = -1;
	}
	Nflag_[columnOut_] = 1;

	// 2. Now we can update the dual
	row_ep.clear();
	row_ap.clear();
	row_ep.count = 1;
	row_ep.index[0] = rowOut_;
	row_ep.array[rowOut_] = 1;
	row_ep.packFlag = true;
	factor.btran(row_ep, row_ep_hist_dsty);
	matrix.price_by_row(row_ap, row_ep);
	row_ep_hist_dsty = 0.95 * row_ep_hist_dsty + 0.05 * row_ep.count / numRow;

	thetaDual_ = dual_[columnIn_] / alpha_;
	for (int i = 0; i < row_ap.count; i++) {
		int iCol = row_ap.index[i];
		dual_[iCol] -= thetaDual_ * row_ap.array[iCol];
	}
	for (int i = 0; i < row_ep.count; i++) {
		int iGet = row_ep.index[i];
		int iCol = iGet + numCol;
		dual_[iCol] -= thetaDual_ * row_ep.array[iGet];
	}

	// Dual for the pivot
	dual_[columnIn_] = 0;
	dual_[columnOut_] = -thetaDual_;

	// Update factor basis
	factor.update(column, row_ep, rowOut_, solving_need_invert);
	matrix.update(columnIn_, columnOut_);
	if (++count_update >= limit_update)
		solving_need_invert = true;
	count_iteration++;

	// Recompute the obj
	objective = 0;
	for (int i = 0; i < numTotal; i++)
		if (Nflag_[i])
			objective += value_[i] * dual_[i];

}

void HModel::sense() {
	/*
	 * Ranging 1.1. set xi and dj strictly feasible
	 */
	vector<double> xi = Bvalue_;
	for (int i = 0; i < numRow; i++) {
		xi[i] = max(xi[i], Blower_[i]);
		xi[i] = min(xi[i], Bupper_[i]);
	}

	vector<double> dj = dual_;
	for (int j = 0; j < numTotal; j++) {
		if (Nflag_[j] && (lower_[j] != upper_[j])) {
			if (value_[j] == lower_[j])
				dj[j] = max(dj[j], 0.0);
			if (value_[j] == upper_[j])
				dj[j] = min(dj[j], 0.0);
			if (lower_[j] == -H_INF && upper_[j] == H_INF)
				dj[j] = 0;
		}
	}

	/*
	 * Ranging 1.2. prepare "delta" space
	 */
	vector<double> dxi_inc(numRow);
	vector<double> dxi_dec(numRow);
	for (int i = 0; i < numRow; i++) {
		dxi_inc[i] = Bupper_[i] - xi[i];
		dxi_dec[i] = Blower_[i] - xi[i];
	}

	vector<double> ddj_inc(numTotal);
	vector<double> ddj_dec(numTotal);
	for (int j = 0; j < numTotal; j++) {
		if (Nflag_[j]) {
			ddj_inc[j] = (value_[j] == lower_[j]) ? +H_INF : -dj[j];
			ddj_dec[j] = (value_[j] == upper_[j]) ? -H_INF : -dj[j];
		}
	}

	/*
	 * Ranging 1.3. prepare "theta" space
	 */
	const double tol_a = 1e-9;
	const double THETA_INF = H_INF / 1e40;

	vector<double> txj_inc(numTotal, +THETA_INF); // theta
	vector<double> axj_inc(numTotal, 0); // alpha
	vector<int> ixj_inc(numTotal, -1); // i-out
	vector<int> wxj_inc(numTotal, 0); // which bound is limiting
	vector<int> jxj_inc(numTotal, -1); // j = n(i), (with bound flip)

	vector<double> txj_dec(numTotal, -THETA_INF);
	vector<double> axj_dec(numTotal, 0);
	vector<int> ixj_dec(numTotal, -1);
	vector<int> wxj_dec(numTotal, 0);
	vector<int> jxj_dec(numTotal, -1);

	vector<double> tci_inc(numRow, +THETA_INF); // theta
	vector<double> aci_inc(numRow, 0); // alpha
	vector<int> jci_inc(numRow, -1); // column index

	vector<double> tci_dec(numRow, -THETA_INF);
	vector<double> aci_dec(numRow, 0);
	vector<int> jci_dec(numRow, -1);

	// Major "theta" loop
	for (int j = 0; j < numTotal; j++) {
		// Skip basic column
		if (!Nflag_[j])
			continue;

		// Form updated column
		column.clear();
		matrix.collect_aj(column, j, 1);
		factor.ftran(column, 0);
		int nWork = 0;
		for (int k = 0; k < column.count; k++) {
			int iRow = column.index[k];
			double alpha = column.array[iRow];
			if (fabs(alpha) > tol_a) {
				iWork_[nWork] = iRow;
				dWork_[nWork] = alpha;
				nWork++;
			}
		}

		// Standard primal ratio test
		double myt_inc = +THETA_INF;
		double myt_dec = -THETA_INF;
		int myk_inc = -1;
		int myk_dec = -1;
		for (int k = 0; k < nWork; k++) {
			int i = iWork_[k];
			double alpha = dWork_[k];
			double theta_inc = (alpha < 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
			double theta_dec = (alpha > 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
			if (myt_inc > theta_inc)
				myt_inc = theta_inc, myk_inc = k;
			if (myt_dec < theta_dec)
				myt_dec = theta_dec, myk_dec = k;
		}

		if (myk_inc != -1) {
			int i = iWork_[myk_inc];
			double alpha = dWork_[myk_inc];
			ixj_inc[j] = i;
			axj_inc[j] = alpha;
			txj_inc[j] = (alpha < 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
			wxj_inc[j] = (alpha < 0 ? +1 : -1);
		}

		if (myk_dec != -1) {
			int i = iWork_[myk_dec];
			double alpha = dWork_[myk_dec];
			ixj_dec[j] = i;
			axj_dec[j] = alpha;
			txj_dec[j] = (alpha > 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
			wxj_dec[j] = (alpha > 0 ? +1 : -1);
		}

		// Accumulated dual ratio test
		double myd_inc = ddj_inc[j];
		double myd_dec = ddj_dec[j];
		for (int k = 0; k < nWork; k++) {
			int i = iWork_[k];
			double alpha = dWork_[k];
			double theta_inc = (alpha < 0 ? myd_inc : myd_dec) / -alpha;
			double theta_dec = (alpha > 0 ? myd_inc : myd_dec) / -alpha;
			if (tci_inc[i] > theta_inc)
				tci_inc[i] = theta_inc, aci_inc[i] = alpha, jci_inc[i] = j;
			if (tci_dec[i] < theta_dec)
				tci_dec[i] = theta_dec, aci_dec[i] = alpha, jci_dec[i] = j;
		}

	}

	// Additional j-out for primal ratio test (considering bound flip)
	for (int j = 0; j < numTotal; j++) {
		if (Nflag_[j]) {
			// J-out for x_j = l_j
			if (Nmove_[j] == +1) {
				double value = value_[j] + txj_inc[j];
				if (ixj_inc[j] != -1 && value <= upper_[j]) {
					jxj_inc[j] = Bindex_[ixj_inc[j]];
				} else if (value > upper_[j]) {
					jxj_inc[j] = j;
				}
			}
			// J-out for x_j = u_j
			if (Nmove_[j] == -1) {
				double value = value_[j] + txj_dec[j];
				if (ixj_dec[j] != -1 && value >= lower_[j]) {
					jxj_dec[j] = Bindex_[ixj_dec[j]];
				} else if (value < lower_[j]) {
					jxj_dec[j] = j;
				}
			}
			// J-out for free variable
			if (lower_[j] == -H_INF && upper_[j] == H_INF) {
				if (ixj_inc[j] != -1)
					jxj_inc[j] = jxj_dec[j] = Bindex_[ixj_inc[j]];
				if (ixj_dec[j] != -1)
					jxj_inc[j] = jxj_dec[j] = Bindex_[ixj_dec[j]];

			}
		}
	}

	/*
	 * Ranging 2. cost ranging
	 */
	vector<double> c_up_c(numTotal), c_dn_c(numTotal);
	vector<double> c_up_f(numTotal), c_dn_f(numTotal);
	vector<int> c_up_e(numTotal), c_dn_e(numTotal);
	vector<int> c_up_l(numTotal), c_dn_l(numTotal);

	/*
	 * Ranging 2.1. non-basic cost ranging
	 */
	for (int j = 0; j < numCol; j++) {
		if (Nflag_[j]) {
			// Primal value and its sign
			double value = value_[j];
			double vsign = (value > 0) ? 1 : (value < 0 ? -1 : 0);

			// Increase c_j
			if (ddj_inc[j] != H_INF) {
				c_up_c[j] = cost_[j] + ddj_inc[j];
				c_up_f[j] = objective + value * ddj_inc[j];
				c_up_e[j] = j;
				c_up_l[j] = jxj_dec[j];
			} else {
				c_up_c[j] = H_INF;
				c_up_f[j] = objective + vsign * H_INF;
				c_up_e[j] = -1;
				c_up_l[j] = -1;
			}

			// Decrease c_j
			if (ddj_dec[j] != H_INF) {
				c_dn_c[j] = cost_[j] + ddj_dec[j];
				c_dn_f[j] = objective + value * ddj_dec[j];
				c_dn_e[j] = j;
				c_dn_l[j] = jxj_inc[j];
			} else {
				c_up_c[j] = -H_INF;
				c_up_f[j] = objective - vsign * H_INF;
				c_up_e[j] = -1;
				c_up_l[j] = -1;
			}
		}
	}

	/*
	 * Ranging 2.2. basic cost ranging
	 */

	for (int i = 0; i < numRow; i++) {
		if (Bindex_[i] < numCol) {
			// Primal variable and its sign
			int j = Bindex_[i], je;
			double value = xi[i];
			double vsign = (value > 0) ? 1 : (value < 0 ? -1 : 0);

			// Increase c_i
			if (jci_inc[i] != -1) {
				c_up_c[j] = cost_[j] + tci_inc[i];
				c_up_f[j] = objective + value * tci_inc[i];
				c_up_e[j] = je = jci_inc[i];
				c_up_l[j] = Nmove_[je] > 0 ? jxj_inc[je] : jxj_dec[je];
			} else {
				c_up_c[j] = H_INF;
				c_up_f[j] = objective + vsign * H_INF;
				c_up_e[j] = -1;
				c_up_l[j] = -1;
			}

			// Decrease c_i
			if (jci_dec[i] != -1) {
				c_dn_c[j] = cost_[j] + tci_dec[i];
				c_dn_f[j] = objective + value * tci_dec[i];
				c_dn_e[j] = je = jci_dec[i];
				c_dn_l[j] = Nmove_[je] > 0 ? jxj_inc[je] : jxj_dec[je];
			} else {
				c_dn_c[j] = -H_INF;
				c_dn_f[j] = objective - H_INF * vsign;
				c_dn_e[j] = -1;
				c_dn_l[j] = -1;
			}
		}
	}

	/*
	 * Ranging 3. bounds ranging
	 */
	vector<double> b_up_b(numTotal), b_dn_b(numTotal);
	vector<double> b_up_f(numTotal), b_dn_f(numTotal);
	vector<int> b_up_e(numTotal), b_dn_e(numTotal);
	vector<int> b_up_l(numTotal), b_dn_l(numTotal);

	/*
	 * Ranging 3.1. non-basic bounds ranging
	 */
	for (int j = 0; j < numTotal; j++) {
		if (Nflag_[j]) {
			// FREE variable
			if (lower_[j] == -H_INF && upper_[j] == H_INF) {
				b_up_b[j] = H_INF;
				b_up_f[j] = objective;
				b_up_e[j] = -1;
				b_up_l[j] = -1;
				b_dn_b[j] = -H_INF;
				b_dn_f[j] = objective;
				b_dn_e[j] = -1;
				b_dn_l[j] = -1;
				continue;
			}

			// Dual value and its sign
			double dualv = dj[j];
			double dsign = (dualv > 0) ? 1 : (dualv < 0 ? -1 : 0);

			// Increase x_j
			if (ixj_inc[j] != -1) {
				int i = ixj_inc[j];
				b_up_b[j] = value_[j] + txj_inc[j];
				b_up_f[j] = objective + txj_inc[j] * dualv;
				b_up_e[j] = wxj_inc[j] > 0 ? jci_inc[i] : jci_dec[i];
				b_up_l[j] = Bindex_[i];
			} else {
				b_up_b[j] = H_INF;
				b_up_f[j] = objective + H_INF * dsign;
				b_up_e[j] = -1;
				b_up_l[j] = -1;
			}

			// Check if b_up_b > upper
			if (value_[j] != upper_[j] && b_up_b[j] > upper_[j]) {
				b_up_b[j] = upper_[j];
				b_up_f[j] = objective + (upper_[j] - lower_[j]) * dualv;
				b_up_e[j] = j;
				b_up_l[j] = j;
			}

			// Decrease x_j
			if (ixj_dec[j] != -1) {
				int i = ixj_dec[j];
				b_dn_b[j] = value_[j] + txj_dec[j];
				b_dn_f[j] = objective + txj_dec[j] * dualv;
				b_dn_e[j] = wxj_dec[j] > 0 ? jci_inc[i] : jci_dec[i];
				b_dn_l[j] = Bindex_[i];
			} else {
				b_dn_b[j] = -H_INF;
				b_dn_f[j] = objective - H_INF * dsign;
				b_dn_e[j] = -1;
				b_dn_l[j] = -1;
			}

			// Check if b_dn_b < lower
			if (value_[j] != lower_[j] && b_dn_b[j] < lower_[j]) {
				b_dn_b[j] = lower_[j];
				b_dn_f[j] = objective + (lower_[j] - upper_[j]) * dualv;
				b_dn_e[j] = j;
				b_dn_l[j] = j;
			}
		}
	}

	/*
	 * Ranging 3.2. basic bounds ranging
	 */
	for (int i = 0; i < numRow; i++) {
		for (int dir = -1; dir <= 1; dir += 2) {
			int j = Bindex_[i];
			double& newx = dir == -1 ? b_dn_b[j] : b_up_b[j];
			double& newf = dir == -1 ? b_dn_f[j] : b_up_f[j];
			int& j_enter = dir == -1 ? b_dn_e[j] : b_up_e[j];
			int& j_leave = dir == -1 ? b_dn_l[j] : b_up_l[j];

			int j_in = dir == -1 ? jci_inc[i] : jci_dec[i];
			double a_in = dir == -1 ? aci_inc[i] : aci_dec[i];
			if (j_in != -1) {
				int jmove = Nmove_[j_in];
				int i_out = jmove > 0 ? ixj_inc[j_in] : ixj_dec[j_in];
				int j_out = jmove > 0 ? jxj_inc[j_in] : jxj_dec[j_in];
				int w_out = jmove > 0 ? wxj_inc[j_in] : wxj_dec[j_in];
				double tt = jmove > 0 ? txj_inc[j_in] : txj_dec[j_in];
				if (j_out == j_in) {
					// Bound flip
					double delta = jmove * (upper_[j_in] - lower_[j_in]);
					newx = xi[i] - delta * a_in;
					newf = objective + delta * dual_[j_in];
					j_enter = j_in;
					j_leave = j_out;
				} else if (j_out != -1) {
					// Regular
					double delta = w_out > 0 ? dxi_inc[i_out] : dxi_dec[i_out];
					double a_out = jmove > 0 ? axj_inc[j_in] : axj_dec[j_in];
					newx = xi[i] + delta * a_in / a_out;
					newf = objective + tt * dual_[j_in];
					j_enter = j_in;
					j_leave = j_out;
				} else {
					// Primal ratio test failed - change unlimitedly
					// While still limited by it's own bounds
					// It's own bounds could just be inf
					newx = dir == -1 ? lower_[j] : upper_[j];
					newf = objective;
					j_enter = -1;
					j_leave = -1;
				}
			} else {
				// Dual ratio test failed - just stay
				newx = xi[i];
				newf = objective;
				j_enter = -1;
				j_leave = -1;
			}
		}
	}

	/*
	 * Ranging 4.1. Scale back
	 */
	for (int j = 0; j < numCol; j++) {
		c_up_c[j] /= (c_up_c[j] == +H_INF) ? 1 : col_scale[j];
		c_dn_c[j] /= (c_dn_c[j] == -H_INF) ? 1 : col_scale[j];
		b_up_b[j] *= (b_up_b[j] == +H_INF) ? 1 : col_scale[j];
		b_dn_b[j] *= (b_dn_b[j] == +H_INF) ? 1 : col_scale[j];

	}
	for (int i = 0, j = numCol; i < numRow; i++, j++) {
		b_up_b[j] /= (b_up_b[j] == +H_INF) ? 1 : row_scale[i];
		b_dn_b[j] /= (b_dn_b[j] == +H_INF) ? 1 : row_scale[i];
	}

	/*
	 * Ranging 4.1.1 Trim small value to zero
	 */
	for (int j = 0; j < numCol; j++) {
		if (fabs(c_up_c[j]) < H_TT)
			c_up_c[j] = 0;
		if (fabs(c_dn_c[j]) < H_TT)
			c_dn_c[j] = 0;
		if (fabs(b_up_b[j]) < H_TT)
			b_up_b[j] = 0;
		if (fabs(b_dn_b[j]) < H_TT)
			b_dn_b[j] = 0;
	}
	for (int i = 0, j = numCol; i < numRow; i++, j++) {
		if (fabs(b_up_b[j]) < H_TT)
			b_up_b[j] = 0;
		if (fabs(b_dn_b[j]) < H_TT)
			b_dn_b[j] = 0;
	}

	/*
	 * Ranging 4.2. Put to output buffer
	 */
	houtput.b_up_b = b_up_b, houtput.b_dn_b = b_dn_b;
	houtput.b_up_f = b_up_f, houtput.b_dn_f = b_dn_f;
	houtput.b_up_e = b_up_e, houtput.b_dn_e = b_dn_e;
	houtput.b_up_l = b_up_l, houtput.b_dn_l = b_dn_l;

	houtput.c_up_c = c_up_c, houtput.c_dn_c = c_dn_c;
	houtput.c_up_f = c_up_f, houtput.c_dn_f = c_dn_f;
	houtput.c_up_e = c_up_e, houtput.c_dn_e = c_dn_e;
	houtput.c_up_l = c_up_l, houtput.c_dn_l = c_dn_l;
}

