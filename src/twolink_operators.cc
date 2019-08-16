/*
 * twolink_operators.cc
 *
 *  Created on: 9 Aug 2019
 *      Author: reisinger
 */

#include <LinkPath.hh>
#include <linear_algebra.hh>
#include <sublattice_algebra.hh>
#include <twolink_operators.hh>

/**
 * function naming:
 * 'x' computes the 'two-link' operator as direct product of the expression to its left with the one to its right
 * 'U' is a single link in temporal direction, at R=0 on the left of 'x' and R=rsep on the right
 * 'C' is the clover, 'low'/'upp' mean only the lower/upper half
 * 'I' is the identity
 */

void U_x_U(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true);

	so_eq_cm_x_cm(result, T0.path, TR.path);
}

void UU_x_UU(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true)(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true)(0, true);

	so_eq_cm_x_cm(result, T0.path, TR.path);
}

void UU_x_UUClow(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
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

	so_eq_cm_x_cm(result, T0.path, U);
}

void UU_x_CuppUU(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
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

	so_eq_cm_x_cm(result, T0.path, U);
}

void UU_x_UCU(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true)(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true);

	double clov[SUN_elems];
	LinkPath plaq(sub_gauge_field, T, L, { t + 1, x, y, z });

	plaq.move(dir, rsep)(0, true)(dir, true)(0, false)(dir, false);
	cm_eq_cm(clov, plaq.path);

	plaq.reset( { t + 1, x, y, z }).move(dir, rsep)(dir, false)(0, true)(dir, true)(0, false);
	cm_pl_eq_cm(clov, plaq.path);

	plaq.reset( { t + 1, x, y, z }).move(dir, rsep)(0, false)(dir, false)(0, true)(dir, true);
	cm_pl_eq_cm(clov, plaq.path);

	plaq.reset( { t + 1, x, y, z }).move(dir, rsep)(dir, true)(0, false)(dir, false)(0, true);
	cm_pl_eq_cm(clov, plaq.path);

	double U[SUN_elems];
	cm_eq_cm_dag(U, clov);
	cm_ti_eq_re(U, -1.0);
	cm_pl_eq_cm(clov, U);
	cm_ti_eq_re(clov, 0.125);

	cm_eq_cm_ti_cm(U, TR.path, clov);
	TR.reset( { t + 1, x, y, z }).move(dir, rsep)(0, true);
	cm_eq_cm_ti_cm(clov, U, TR.path);

	so_eq_cm_x_cm(result, T0.path, clov);
}

void IU_x_IU(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t + 1, x, y, z });
	T0(0, true);
	LinkPath TR(sub_gauge_field, T, L, { t + 1, x, y, z });
	TR.move(dir, rsep)(0, true);
	so_eq_cm_x_cm(result, T0.path, TR.path);
}

void UI_x_UI(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true);
	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true);
	so_eq_cm_x_cm(result, T0.path, TR.path);
}
