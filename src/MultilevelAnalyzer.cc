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

#include <fields.hh>
#include <io.hh>
#include <ranlux.hh>
#include <heatbath.hh>

#include <global_defs.hh>

#include <sublattice_algebra.hh>
#include <sublattice_fields.hh>

#include <MultilevelAnalyzer.hh>

//public

MultilevelAnalyzer::MultilevelAnalyzer(
		int T, int L, int WL_R, std::vector<int> level_thickness,
		std::string config_prefix, std::vector<int> level_config_num,
		std::vector<void (*)(double*, const double*, int, int, int&, int, int, int, int, int)> lowest_level_functions,
		std::vector<std::vector<std::vector<int> > > field_compositions,
		bool generate_configs, double beta, int seed, std::vector<int> level_updates, bool save) :
		T(T), L(L), WL_R(WL_R), level_thickness(level_thickness),
				config_prefix(config_prefix), level_config_num(level_config_num),
				field_compositions(field_compositions), lowest_level_functions(lowest_level_functions),
				generate_configs(generate_configs), beta(beta), seed(seed), level_updates(level_updates), save(save)
{
	if (generate_configs) {
		InitializeRand(seed);
		Gauge_Field_Alloc_silent(&config_buf, T, L);
	}
}

MultilevelAnalyzer::~MultilevelAnalyzer() {
	if (generate_configs)
		Gauge_Field_Free(&config_buf);
}

std::string MultilevelAnalyzer::tag_to_string(const std::vector<int>& tag) {
	std::ostringstream tag_oss;
	for (int i = 0; i < tag.size(); ++i)
		tag_oss << "." << std::setfill('0') << std::setw(log10(level_config_num.at(i)) + 1) << tag.at(i);
	return tag_oss.str();
}

std::string MultilevelAnalyzer::config_filename(const std::vector<int>& tag) {
	std::ostringstream filename_oss;
	filename_oss << config_prefix << tag_to_string(tag);
	return filename_oss.str();
}

void MultilevelAnalyzer::compute_sublattice_fields(const std::vector<int>& conf_tag, const int level, double** T_fields) {
	const bool is_lowest = (level == field_compositions.size() - 1);

	const int config_num = (level == 0 ? 1 : level_config_num[level]);
	const int timeslice_thickness = level_thickness[level == 0 ? 1 : level];
	const int lower_level_field_num = (is_lowest ? 1 : field_compositions[level + 1].size());

	for (int conf = 1; conf <= config_num; ++conf) {
		std::vector<int> curr_tag(conf_tag);
		if (level != 0)
			curr_tag.push_back(conf);

		if (generate_configs) {
			if (level == 0)
				read_gauge_field(config_buf, config_filename(curr_tag).c_str(), T, L);
			update_sublattice_gauge_field(curr_tag);
			if (level != 0 && save)
				write_config(curr_tag);
		}

		double* lower_level_fields[lower_level_field_num];
		if (is_lowest)
			obtain_sublattice_gauge_field(lower_level_fields[0], curr_tag);
		else {
			std::cerr << "Allocating sublattice fields on config " << tag_to_string(curr_tag) << " ... ";
			for (int i = 0; i < lower_level_field_num; ++i)
				T_field_alloc_zero(lower_level_fields[i], 3, T / level_thickness[level + 1], L);
			std::cerr << "ok\n";
			compute_sublattice_fields(curr_tag, level + 1, lower_level_fields);
		}

		for (int t = 0; t < T; t += timeslice_thickness) {
			for (int x = 0; x < L; ++x) {
				for (int y = 0; y < L; ++y) {
					for (int z = 0; z < L; ++z) {
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
												operator_component, lower_level_fields[0], T, L, curr_t, x, y, z, i, WL_R);
									} else {
										operator_component = lower_level_fields[lower_level_field_index]
												+ T_field_index(curr_t, x, y, z, 3, i - 1, T, L, level_thickness[level + 1]);
									}

									double Ttemp[SO_elems];

									so_eq_so_ti_so(Ttemp, curr_operator, operator_component);

									if (is_lowest)
										delete[] operator_component;
									so_eq_so(curr_operator, Ttemp);

									if (!is_lowest) //at lowest level, curr_t is incremented by the twolink_operators function
										curr_t += level_thickness[level + 1];
								}

								so_pl_eq_so(T_fields[curr_field_index] + T_field_index(t, x, y, z, 3, i - 1, T, L, timeslice_thickness),
										curr_operator);
							}

						}
					}
				}
			}
		}
		for (int i = 0; i < lower_level_field_num; ++i)
			level_field_free(lower_level_fields[i], level);
	}
	for (int i = 0; i < field_compositions[level].size(); ++i)
		T_field_di_eq_re(T_fields[i], config_num, 3, T, L, timeslice_thickness);
}

//private

void MultilevelAnalyzer::obtain_sublattice_gauge_field(double*& sub_gauge_field, const std::vector<int>& tag) {
	if (generate_configs)
		sub_gauge_field = config_buf;
	else {
		std::cerr << "Reading config " << tag_to_string(tag) << " ... ";
		Gauge_Field_Alloc_silent(&sub_gauge_field, T, L);
		read_gauge_field(sub_gauge_field, config_filename(tag).c_str(), T, L);
		std::cerr << "ok\n";
	}
}

void MultilevelAnalyzer::update_sublattice_gauge_field(const std::vector<int>& tag) {
	std::cerr << "Generating config " << tag_to_string(tag) << " ... ";
	std::vector<int> boundary_ts;
	int level = tag.size() - 1;
	if (level != 0)
		for (int t = 0; t < T; t += level_thickness[level])
			boundary_ts.push_back(t);
	for (int i_swp = 0; i_swp < level_updates[level]; ++i_swp)
		do_sweep(config_buf, T, L, beta, boundary_ts);
	std::cerr << "ok\n";
}

void MultilevelAnalyzer::write_config(const std::vector<int>& tag) {
	std::ostringstream config_filename_oss;
	config_filename_oss << config_prefix << ".multilevel" << tag_to_string(tag);
	std::ostringstream header_oss;
	header_oss << "generated during multilevel : " << beta << " " << T << " " << L << " " << seed;
	write_gauge_field(config_buf, config_filename_oss.str().c_str(), T, L, header_oss.str().c_str());
}

void MultilevelAnalyzer::level_field_free(double*& level_field, int level) {
	if (level == field_compositions.size() - 1) {
		if (!generate_configs)
			Gauge_Field_Free(&level_field);
	} else
		T_field_free(level_field);
}
