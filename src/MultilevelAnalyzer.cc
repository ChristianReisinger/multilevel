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

//private

void MultilevelAnalyzer::obtain_sublattice_gauge_field(double*& sub_gauge_field, const std::vector<int>& tag) {
	Gauge_Field_Alloc(&sub_gauge_field, T, L);
	read_gauge_field(sub_gauge_field, config_filename(tag).c_str(), T, L);
}
