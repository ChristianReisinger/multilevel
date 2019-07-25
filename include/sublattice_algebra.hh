/*
 * sublattice_algebra.hh
 *
 *  Created on: 24 Jul 2019
 *      Author: reisinger
 */

#include <global_defs.hh>
#include <linear_algebra.hh>

#ifndef INCLUDE_SUBLATTICE_ALGEBRA_HH_
#define INCLUDE_SUBLATTICE_ALGEBRA_HH_

constexpr int SO_elems = 2 * SUN_N * SUN_N * SUN_N * SUN_N;

int so_superindex(int alpha, int beta, int gamma, int delta) {
	return 2 * (delta + SUN_elems * (gamma + SUN_elems * (beta + SUN_elems * alpha)));
}

int cm_superindex(int row, int col) {
	return 2 * (SUN_N * row + col);
}

void so_eq_zero(double* const result) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] = 0.0;
}

void so_eq_id(double* const result) {
	so_eq_zero(result);
	for (int a = 0; a < SUN_N; ++a) {
		for (int b = 0; b < SUN_N; ++b) {
			result[so_superindex(a, a, b, b)] = 1;
		}
	}
}

void so_eq_so(double* const result, const double* const T) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] = T[i];
}

void so_eq_so_ti_so(double* const result, const double* const T1, const double* const T2) {
	for (int a = 0; a < SUN_N; ++a) {
		for (int b = 0; b < SUN_N; ++b) {
			for (int c = 0; c < SUN_N; ++c) {
				for (int d = 0; d < SUN_N; ++d) {
					double* result_elem_ptr = result + so_superindex(a, b, c, d);
					for (int i = 0; i < SUN_N; ++i) {
						for (int j = 0; j < SUN_N; ++j) {
							double* T1_elem_ptr = T1 + so_superindex(a, i, c, j);
							double* T2_elem_ptr = T2 + so_superindex(i, b, j, d);
							complex T1_elem { *T1_elem_ptr, *(T1_elem_ptr + 1) };
							complex T2_elem { *T2_elem_ptr, *(T2_elem_ptr + 1) };
							complex result_elem;
							co_eq_co_ti_co(&result_elem, &T1_elem, &T2_elem);
							*result_elem_ptr += result_elem.re;
							*(result_elem_ptr + 1) += result_elem.im;
						}
					}
				}
			}
		}
	}
}

void so_pl_eq_so(double* const result, const double* const T) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] += T[i];
}

void so_ti_eq_re(double* const result, double d) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] *= d;
}

/**
 * compute the sublattice operator from the "two link" operators T0 at r = 0 and TR at r = R = spatial Wilson loop size
 */
void so_eq_cm_x_cm(double* const result, const double* const T0, const double* const TR) {
	double T0_dag[SUN_elems];
	cm_eq_cm_dag(T0_dag, T0);
	for (int a = 0; a < SUN_N; ++a) {
		for (int b = 0; b < SUN_N; ++b) {
			for (int c = 0; c < SUN_N; ++c) {
				for (int d = 0; d < SUN_N; ++d) {
					double* result_elem_ptr = result + so_superindex(a, b, c, d);
					double* T0_elem_ptr = T0_dag + cm_superindex(b, a);
					double* TR_elem_ptr = TR + cm_superindex(c, d);
					complex T0_elem { *T0_elem_ptr, *(T0_elem_ptr + 1) };
					complex TR_elem { *TR_elem_ptr, *(TR_elem_ptr + 1) };
					complex result_elem;
					co_eq_co_ti_co(&result_elem, &T0_elem, &result_elem);
					*result_elem_ptr = result_elem.re;
					*(result_elem_ptr + 1) = result_elem.im;
				}
			}
		}
	}
}

#endif /* INCLUDE_SUBLATTICE_ALGEBRA_HH_ */
