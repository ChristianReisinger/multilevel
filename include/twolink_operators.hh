/*
 * twolink_operators.hh
 *
 *  Created on: 6 Aug 2019
 *      Author: reisinger
 */

#ifndef INCLUDE_TWOLINK_OPERATORS_HH_
#define INCLUDE_TWOLINK_OPERATORS_HH_

void U_x_U(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void UU_x_UU(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void UU_x_UUClow(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void UU_x_CuppUU(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void UU_x_UCU(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void IU_x_IU(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void UI_x_UI(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

#endif /* INCLUDE_TWOLINK_OPERATORS_HH_ */
