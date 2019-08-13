/*
 * MultilevelAnalyzer.cc
 *
 *  Created on: 9 Aug 2019
 *      Author: reisinger
 */

#include <string>
#include <vector>
#include <iomanip>
#include <sstream>

#include <fields.hh>
#include <io.hh>
#include <sublattice_algebra.hh>
#include <sublattice_fields.hh>
#include <global_defs.hh>

#include <MultilevelAnalyzer.hh>

//public

MultilevelAnalyzer::MultilevelAnalyzer(
		int T, int L, int WL_R, std::vector<int> level_thickness,
		std::string config_prefix, std::vector<int> level_config_num,
		std::vector<void (*)(double*, const double*, int, int, int, int, int, int, int, int)> lowest_level_functions,
		std::vector<std::vector<std::vector<int> > > field_compositions) :
		T(T), L(L), WL_R(WL_R), level_thickness(level_thickness),
				config_prefix(config_prefix), level_config_num(level_config_num),
				field_compositions(field_compositions), lowest_level_functions(lowest_level_functions) {
}

std::string MultilevelAnalyzer::config_filename(const std::vector<int>& tag) {
	std::ostringstream filename_oss;
	filename_oss << config_prefix;
	for (int i = 0; i < tag.size(); ++i)
		filename_oss << "." << std::setfill('0') << std::setw(log10(level_config_num.at(i)) + 1) << tag.at(i);
	return filename_oss.str();
}

/**
 * @param field_compositions a vector V at index L in field_compositions describes compositions
 * 			of the fields F_i(L) at level L from those at level L+1.
 *			A vector W at index i in V contains indices j of the fields at level L+1.
 *			The field F_i(L) is obtained by multiplying fields F_j(L+1) for all j in the order
 *			they appear in W.
 *			At the lowest level (= field_compositions.size() - 1), j instead is an index of a
 *			function in lowest_level_functions.
 */
void MultilevelAnalyzer::compute_sublattice_fields(const std::vector<int>& conf_tag, int level, double** T_fields) {
	const bool is_lowest = (level == field_compositions.size() - 1);

	const int config_num = (level == 0 ? 1 : level_config_num[level]);
	const int timeslice_thickness = level_thickness[level + (level == 0 ? 1 : 0)];
	const int lower_level_field_num = (is_lowest ? 1 : field_compositions[level + 1].size());

	for (int conf = 1; conf <= config_num; ++conf) {
		std::vector<int> curr_tag(conf_tag);
		if (level != 0)
			curr_tag.push_back(conf);

		double* lower_level_fields[lower_level_field_num];
		if (is_lowest)
			obtain_sublattice_gauge_field(lower_level_fields[0], curr_tag);
		else {
			for (int i = 0; i < lower_level_field_num; ++i)
				T_field_alloc_zero(lower_level_fields[i], 3, T / level_thickness[level + 1], L);
			compute_sublattice_fields(curr_tag, level + 1, lower_level_fields);
		}

		for (int t = 0; t < T; t += timeslice_thickness) {
			for (int x = 0; x < L; ++x) {
				for (int y = 0; y < L; ++y) {
					for (int z = 0; z < L; ++z) {
						for (int i = 0; i < 3; ++i) {

							for (int curr_field_index = 0; curr_field_index < field_compositions[level].size(); ++curr_field_index) {
								double curr_operator[SO_elems];
								so_eq_id(curr_operator);
								int curr_t = t;

								for (int lower_level_field_index : field_compositions[level][curr_field_index]) {
									double* component_operator;
									if (is_lowest) {
										component_operator = new double[SO_elems];
										lowest_level_functions[curr_field_index](
												component_operator, lower_level_fields[0], T, L, t, x, y, z, i + 1, WL_R);
									} else {
										component_operator = lower_level_fields[lower_level_field_index]
												+ T_field_index(curr_t, x, y, z, 3, i, T, L, level_thickness[level + 1]);
									}

									double Ttemp[SO_elems];
									so_eq_so_ti_so(Ttemp, curr_operator, component_operator);
									if (is_lowest)
										delete[] component_operator;
									so_eq_so(curr_operator, Ttemp);

									curr_t += level_thickness[level + 1];
								}

								so_pl_eq_so(T_fields[curr_field_index] + T_field_index(t, x, y, z, 3, i, T, L, timeslice_thickness),
										curr_operator);
							}

						}
					}
				}
			}
		}
		for (int i = 0; i < lower_level_field_num; ++i)
			if (is_lowest)
				Gauge_Field_Free(&lower_level_fields[i]);
			else
				delete[] lower_level_fields[i];
	}
	for (int i = 0; i < field_compositions[level].size(); ++i)
		T_field_di_eq_re(T_fields[i], config_num, 3, T, L, timeslice_thickness);
}

//private

void MultilevelAnalyzer::obtain_sublattice_gauge_field(double*& sub_gauge_field, const std::vector<int>& tag) {
	Gauge_Field_Alloc(&sub_gauge_field, T, L);
	read_gauge_field(sub_gauge_field, config_filename(tag).c_str(), T, L);
}
