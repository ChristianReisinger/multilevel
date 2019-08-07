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

inline int so_superindex(int alpha, int beta, int gamma, int delta) {
	return 2 * (delta + SUN_N * (gamma + SUN_N * (beta + SUN_N * alpha)));
}

inline int cm_superindex(int row, int col) {
	return 2 * (SUN_N * row + col);
}

void so_ti_eq_re(double* result, double d) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] *= d;
}

void so_eq_zero(double* result) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] = 0.0;
}

void so_eq_id(double* result) {
	so_eq_zero(result);
	for (int a = 0; a < SUN_N; ++a)
		for (int b = 0; b < SUN_N; ++b)
			result[so_superindex(a, a, b, b)] = 1;
}

void so_eq_so(double* result, const double* T) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] = T[i];
}

void so_pl_eq_so(double* result, const double* T) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] += T[i];
}

void so_eq_so_ti_so(double* result, const double* T1, const double* T2) {
	for (int a = 0; a < SUN_N; ++a) {
		for (int b = 0; b < SUN_N; ++b) {
			for (int c = 0; c < SUN_N; ++c) {
				for (int d = 0; d < SUN_N; ++d) {
					double* result_re = result + so_superindex(a, b, c, d);
					double* result_im = result_re + 1;
					*result_re = 0.0;
					*result_im = 0.0;

					for (int i = 0; i < SUN_N; ++i) {
						for (int j = 0; j < SUN_N; ++j) {
							const double* T1_re = T1 + so_superindex(a, i, c, j);
							const double* T1_im = T1_re + 1;

							const double* T2_re = T2 + so_superindex(i, b, j, d);
							const double* T2_im = T2_re + 1;

							*result_re += *T1_re * *T2_re - *T1_im * *T2_im;
							*result_im += *T1_re * *T2_im + *T1_im * *T2_re;
						}
					}
				}
			}
		}
	}
}

/**
 * compute the sublattice operator from the "two link" operators T0 at r = 0 and TR at r = R = spatial Wilson loop size
 */
void so_eq_cm_x_cm(double* result, const double* T0, const double* TR) {
	for (int a = 0; a < SUN_N; ++a) {
		for (int b = 0; b < SUN_N; ++b) {
			for (int c = 0; c < SUN_N; ++c) {
				for (int d = 0; d < SUN_N; ++d) {
					double* result_re = result + so_superindex(a, b, c, d);
					double* result_im = result_re + 1;

					const double* T0_re = T0 + cm_superindex(a, b);
					const double* T0_im = T0_re + 1;

					const double* TR_re = TR + cm_superindex(c, d);
					const double* TR_im = TR_re + 1;

					*result_re = *T0_re * *TR_re + *T0_im * *TR_im;
					*result_im = *T0_re * *TR_im - *T0_im * *TR_re;
				}
			}
		}
	}
}

/**
 * compute the Wilson loop from the sublattice operator SO and spatial Wilson lines
 * S0 at t = 0, and ST at t = T = temporal Wilson loop size
 */
void close_Wilson_loop(complex* WL, const double* SO, const double* S0, const double* ST) {
	co_eq_zero(WL);

	for (int a = 0; a < SUN_N; ++a) {
		for (int b = 0; b < SUN_N; ++b) {
			for (int c = 0; c < SUN_N; ++c) {
				for (int d = 0; d < SUN_N; ++d) {
					const double* S0_re = S0 + cm_superindex(a, c);
					const double* S0_im = S0_re + 1;

					const double* ST_re = ST + cm_superindex(b, d);
					const double* ST_im = ST_re + 1;

					const double* SO_re = SO + so_superindex(a, b, c, d);
					const double* SO_im = SO_re + 1;

					WL->re += *S0_re * *ST_re * *SO_re
							+ *S0_re * *ST_im * *SO_im
							- *S0_im * *ST_re * *SO_im
							+ *S0_im * *ST_im * *SO_re;

					WL->im += *S0_im * *ST_im * *SO_im
							+ *S0_im * *ST_re * *SO_re
							- *S0_re * *ST_im * *SO_re
							+ *S0_re * *ST_re * *SO_im;
				}
			}
		}
	}

	WL->re /= SUN_N;
	WL->im /= SUN_N;
}

#endif /* INCLUDE_SUBLATTICE_ALGEBRA_HH_ */
