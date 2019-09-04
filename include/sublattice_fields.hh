/*
 * sublattice_fields.hh
 *
 *  Created on: 9 Aug 2019
 *      Author: reisinger
 */

#include <geometry.hh>
#include <sublattice_algebra.hh>

#ifndef INCLUDE_SUBLATTICE_FIELDS_HH_
#define INCLUDE_SUBLATTICE_FIELDS_HH_

inline void T_field_alloc_zero(double*& T_field, int n, int timeslice_num, int L) {
	T_field = new double[n * SO_elems * L * L * L * timeslice_num]();
}

inline void T_field_free(double*& T_field) {
	delete[] T_field;
}

inline unsigned long int T_field_index(int t, int x, int y, int z, int n, int i, int T, int L, int timeslice_thickness) {
	return (n * get_index(t / timeslice_thickness, x, y, z, T / timeslice_thickness, L) + i) * SO_elems;
}

inline void T_field_di_eq_re(double* T_field, double re, int n, int T, int L, int timeslice_thickness) {
	for (int i = 0; i < n * SO_elems * L * L * L * T / timeslice_thickness; ++i)
		T_field[i] /= re;
}

#endif /* INCLUDE_SUBLATTICE_FIELDS_HH_ */
