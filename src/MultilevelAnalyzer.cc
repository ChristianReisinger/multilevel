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
#include <iostream>
#include <string>
#include <chrono>

#include <fields.hh>
#include <io.hh>
#include <ranlux.hh>
#include <heatbath.hh>

#include <global_defs.hh>

#include <sublattice_algebra.hh>
#include <sublattice_fields.hh>

#include <MultilevelConfig.hh>
#include <MultilevelAnalyzer.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

//public

MultilevelAnalyzer::MultilevelAnalyzer(MultilevelConfig& multilevel_config, std::set<int> WL_Rs,
		std::vector<std::vector<std::vector<int> > > field_compositions,
		std::vector<void (*)(double*, const double*, int, int, int&, int, int, int, int, int)> lowest_level_functions) :
		config(&multilevel_config), WL_Rs(WL_Rs), field_compositions(field_compositions), lowest_level_functions(lowest_level_functions),
		time_spent_computing_operators(0) {
}

void MultilevelAnalyzer::compute_sublattice_fields(std::map<int, double**> T_fields, const int level) {
	if (top_level && level != 0)
		throw std::invalid_argument("multilevel top level must be 0");
	else
		top_level = false;

	const bool is_lowest = (level == field_compositions.size() - 1);

	const int config_num = (level == 0 ? 1 : config->config_num(level));
	const int timeslice_thickness = config->thickness(level == 0 ? 1 : level);
	const int lower_level_field_num = (is_lowest ? 1 : field_compositions[level + 1].size());

	for (int conf = 1; conf <= config_num; ++conf) {
		config->update(level);

		std::map<int, double**> lower_level_fields;
		double* lowest_level_config;
		if (is_lowest)
			config->get(lowest_level_config);
		else {
			std::cerr << "Allocating sublattice fields on config '" << config->config_filename() << "' ... ";
			for (const auto WL_R : WL_Rs) {
				lower_level_fields[WL_R] = new double*[lower_level_field_num];
				for (int i = 0; i < lower_level_field_num; ++i)
					T_field_alloc_zero(lower_level_fields.at(WL_R)[i], 3, config->T / config->thickness(level + 1), config->L);
			}
			std::cerr << "ok\n";
			compute_sublattice_fields(lower_level_fields, level + 1);
		}

		auto start_time = std::chrono::steady_clock::now();
		for (const auto WL_R : WL_Rs) {
			for (int t = 0; t < config->T; t += timeslice_thickness) {
				for (int x = 0; x < config->L; ++x) {
					for (int y = 0; y < config->L; ++y) {
						for (int z = 0; z < config->L; ++z) {
							for (int i = 1; i < 4; ++i) {

								for (int curr_field_index = 0; curr_field_index < field_compositions[level].size(); ++curr_field_index) {
									double curr_operator[SO_elems];
									so_eq_id(curr_operator);
									int curr_t = t;

									for (int lower_level_field_index : field_compositions[level][curr_field_index]) {
										double* operator_component;
										if (is_lowest) {
											operator_component = new double[SO_elems];
											lowest_level_functions[lower_level_field_index](
													operator_component, lowest_level_config, config->T, config->L,
													curr_t, x, y, z, i, WL_R);
										} else {
											operator_component = lower_level_fields.at(WL_R)[lower_level_field_index]
													+ T_field_index(curr_t, x, y, z, 3, i - 1, config->T, config->L,
															config->thickness(level + 1));
										}

										double Ttemp[SO_elems];

										so_eq_so_ti_so(Ttemp, curr_operator, operator_component);

										if (is_lowest)
											delete[] operator_component;
										so_eq_so(curr_operator, Ttemp);

										if (!is_lowest) //at lowest level, curr_t is incremented by the twolink_operators function
											curr_t += config->thickness(level + 1);
									}

									so_pl_eq_so(T_fields.at(WL_R)[curr_field_index]
											+ T_field_index(t, x, y, z, 3, i - 1, config->T, config->L, timeslice_thickness),
											curr_operator);
								}

							}
						}
					}
				}
			}

			if (!is_lowest)
				for (int i = 0; i < lower_level_field_num; ++i)
					T_field_free(lower_level_fields.at(WL_R)[i]);
		}
		time_spent_computing_operators += std::chrono::steady_clock::now() - start_time;

		if (!is_lowest)
			for (auto& e : lower_level_fields)
				delete[] e.second;
	}
	for (const auto WL_R : WL_Rs)
		for (int i = 0; i < field_compositions[level].size(); ++i)
			T_field_di_eq_re(T_fields.at(WL_R)[i], config_num, 3, config->T, config->L, timeslice_thickness);

	if (level == 0)
		top_level = true;
}

int MultilevelAnalyzer::milliseconds_spent_computing() {
	return std::chrono::duration_cast<std::chrono::milliseconds>(time_spent_computing_operators).count();
}

}
}
}
