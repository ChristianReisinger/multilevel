/*
 * twolink_operators.hh
 *
 *  Created on: 6 Aug 2019
 *      Author: reisinger
 */

#include <LinkPath.hh>
#include <sublattice_algebra.hh>

extern int T;
extern int L;

#ifndef INCLUDE_TWOLINK_OPERATORS_HH_
#define INCLUDE_TWOLINK_OPERATORS_HH_

void compute_T_ti_T(double* T, double* sub_gauge_field, int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true)(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true)(0, true);

	so_eq_cm_x_cm(T, T0.path, TR.path);
}

void compute_T_ti_Tclov_lower_half(double* T, double* sub_gauge_field, int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true)(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true)(0, true);

	double clov[SUN_elems];
	LinkPath plaq(sub_gauge_field, T, L, { t + 2, x, y, z });

	plaq.move(dir, rsep)(dir, true)(0, false)(dir, false)(0, true);
	cm_eq_cm(clov, plaq.path);

	plaq.reset( { t + 2, x, y, z }).move(dir, rsep)(0, false)(dir, false)(0, true)(dir, true);
	cm_pl_eq_cm(clov, plaq.path);

	double U[SUN_elems];
	cm_eq_cm_dag(U, clov);
	cm_ti_eq_re(U, -1.0);
	cm_pl_eq_cm(clov, U);
	cm_ti_eq_re(clov, 0.5);

	cm_eq_cm_ti_cm(U, TR.path, clov);

	so_eq_cm_x_cm(T, T0.path, U);
}

void compute_Tclov_upper_half_ti_T(double* T, double* sub_gauge_field, int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true)(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true)(0, true);

	double clov[SUN_elems];
	LinkPath plaq(sub_gauge_field, T, L, { t, x, y, z });

	plaq.move(dir, rsep)(0, true)(dir, true)(0, false)(dir, false);
	cm_eq_cm(clov, plaq.path);

	plaq.reset( { t, x, y, z }).move(dir, rsep)(dir, false)(0, true)(dir, true)(0, false);
	cm_pl_eq_cm(clov, plaq.path);

	double U[SUN_elems];
	cm_eq_cm_dag(U, clov);
	cm_ti_eq_re(U, -1.0);
	cm_pl_eq_cm(clov, U);
	cm_ti_eq_re(clov, 0.5);

	cm_eq_cm_ti_cm(U, clov, TR.path);

	so_eq_cm_x_cm(T, T0.path, U);
}

#endif /* INCLUDE_TWOLINK_OPERATORS_HH_ */
