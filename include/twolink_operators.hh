/*
 * twolink_operators.hh
 *
 *  Created on: 6 Aug 2019
 *      Author: reisinger
 */

#ifndef INCLUDE_TWOLINK_OPERATORS_HH_
#define INCLUDE_TWOLINK_OPERATORS_HH_

void compute_T_ti_T(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep);

void compute_T_ti_Tclov_lower_half(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep);

void compute_Tclov_upper_half_ti_T(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep);

#endif /* INCLUDE_TWOLINK_OPERATORS_HH_ */
