#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>
#include <chrono>
#include <stdexcept>

#include <fields.hh>
#include <io.hh>
#include <ranlux.hh>
#include <heatbath.hh>

#include <global_defs.hh>

#include <sublattice_algebra.hh>
#include <T_field.hh>

#include <MultilevelConfig.hh>
#include <MultilevelAnalyzer.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

//public

MultilevelAnalyzer::MultilevelAnalyzer(MultilevelConfig& multilevel_config, std::set<int> WL_Rs,
		std::vector<std::map<std::string, std::vector<std::string> > > level_operator_factors,
		std::vector<std::map<std::string, std::vector<bool> > > level_operator_timeslice_defined,
		std::map<std::string, void (*)(double*, const double*, int, int, int&, int, int, int, int, int)> lowest_level_functions) :
		config(&multilevel_config), WL_Rs(WL_Rs),
				level_operator_factors(level_operator_factors), level_operator_timeslice_defined(level_operator_timeslice_defined),
				lowest_level_functions(lowest_level_functions),
				time_spent_computing_operators(0) {
}

std::map<std::string, std::map<int, T_field> > MultilevelAnalyzer::compute_T_fields() {
	std::map<std::string, std::map<int, T_field> > T_fields;
	alloc_T_fields(T_fields, 0);
	compute_sublattice_fields(T_fields, 0);
	config->update(0);
	return T_fields;
}

int MultilevelAnalyzer::milliseconds_spent_computing() {
	return std::chrono::duration_cast<std::chrono::milliseconds>(time_spent_computing_operators).count();
}

//private

void MultilevelAnalyzer::compute_sublattice_fields(std::map<std::string, std::map<int, T_field> >& T_fields, const int level) {
	const bool is_lowest = (level == level_operator_factors.size() - 1);
	const int config_num = (level == 0 ? 1 : config->config_num(level));

	for (int conf = 1; conf <= config_num; ++conf) {
		config->update(level);

		std::map<std::string, std::map<int, T_field> > lower_level_fields;
		double* lowest_level_gauge_field;
		if (is_lowest)
			config->get(lowest_level_gauge_field);
		else {
			std::cerr << "Allocating sublattice fields on config '" << config->config_filename() << "' ... ";
			alloc_T_fields(lower_level_fields, level + 1);
			std::cerr << "ok\n";
			compute_sublattice_fields(lower_level_fields, level + 1);
		}

		auto start_time = std::chrono::steady_clock::now();
		for (auto& name_rfields : T_fields) {
			for (auto& r_fields : name_rfields.second) {
				int WL_R = r_fields.first;
				try {
					for (int t : name_rfields.second.at(WL_R).defined_ts()) {
						for (int x = 0; x < config->L; ++x) {
							for (int y = 0; y < config->L; ++y) {
								for (int z = 0; z < config->L; ++z) {
									for (int i = 1; i < 4; ++i) {
										double curr_operator[SO_elems];
										so_eq_id(curr_operator);
										int curr_t = t;

										for (const std::string& factor_name : level_operator_factors.at(level).at(name_rfields.first)) {
											double* factor;
											if (is_lowest) {
												factor = new double[SO_elems];
												lowest_level_functions.at(factor_name)(
														factor, lowest_level_gauge_field, config->T, config->L, curr_t, x, y, z, i, WL_R);
												//curr_t is incremented by lowest_level_functions
											} else {
												factor = lower_level_fields.at(factor_name).at(WL_R).T_at(curr_t, x, y, z, i);
												curr_t += lower_level_fields.at(factor_name).at(WL_R).timeslice_thickness();
											}

											double Ttemp[SO_elems];
											so_eq_so_ti_so(Ttemp, curr_operator, factor);
											if (is_lowest)
												delete[] factor;

											so_eq_so(curr_operator, Ttemp);
										}

										so_pl_eq_so(name_rfields.second.at(WL_R).T_at(t, x, y, z, i), curr_operator);
									}

								}
							}
						}
					}
				} catch (std::out_of_range& e) {
					throw std::runtime_error("invalid definition of operator '" + name_rfields.first + "'");
				}

			}
		}
		time_spent_computing_operators += std::chrono::steady_clock::now() - start_time;
	}
	for (auto& name_rfields : T_fields)
		for (auto& r_field : name_rfields.second)
			r_field.second /= config_num;
}

void MultilevelAnalyzer::alloc_T_fields(std::map<std::string, std::map<int, T_field> >& T_fields, const int level) {
	for (const auto& fieldname_factors : level_operator_factors.at(level)) {
		const std::string& name = fieldname_factors.first;
		const std::vector<int> level_thickness = config->thickness(level == 0 ? 1 : level);

		std::vector<std::pair<int, bool> > thickness_defined;
		for (int i = 0; i < level_thickness.size(); ++i)
			thickness_defined.push_back( { level_thickness.at(i),
					level_operator_timeslice_defined.at(level).at(name).at(i) });

		std::map<int, T_field> r_fields;
		for (int WL_R : WL_Rs)
			r_fields.insert( { WL_R, T_field(thickness_defined, config->T, config->L) });

		T_fields.insert( { name, r_fields });
	}
}

}
}
}
