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

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

using latticetools_0719::SUN_N;

constexpr int SO_elems = 2 * SUN_N * SUN_N * SUN_N * SUN_N;

inline int so_superindex(int alpha, int beta, int gamma, int delta) {
	return 2 * (delta + SUN_N * (gamma + SUN_N * (beta + SUN_N * alpha)));
}

inline int cm_superindex(int row, int col) {
	return 2 * (SUN_N * row + col);
}

inline void so_ti_eq_re(double* result, double d) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] *= d;
}

inline void so_eq_zero(double* result) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] = 0.0;
}

inline void so_eq_id(double* result) {
	so_eq_zero(result);
	for (int a = 0; a < SUN_N; ++a)
		for (int b = 0; b < SUN_N; ++b)
			result[so_superindex(a, a, b, b)] = 1;
}

inline void so_eq_so(double* result, const double* T) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] = T[i];
}

inline void so_pl_eq_so(double* result, const double* T) {
	for (int i = 0; i < SO_elems; ++i)
		result[i] += T[i];
}

void so_eq_so_ti_so(double* result, const double* T1, const double* T2);

/**
 * compute the sublattice operator from the (path of) temporal links T0 at r = 0 and TR at r = R = spatial Wilson loop size
 */
void so_eq_cm_x_cm(double* result, const double* T0, const double* TR);

/**
 * compute the Wilson loop from the sublattice operator SO and spatial Wilson lines
 * S0 at t = 0, and ST at t = T = temporal Wilson loop size
 */
void close_Wilson_loop(complex* WL, const double* SO, const double* S0, const double* ST);

}
}
}

#endif /* INCLUDE_SUBLATTICE_ALGEBRA_HH_ */
