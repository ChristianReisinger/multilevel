/*
 * MultilevelAnalyzer.hh
 *
 *  Created on: 9 Aug 2019
 *      Author: reisinger
 */

#include <vector>
#include <string>

#ifndef INCLUDE_MULTILEVELANALYZER_HH_
#define INCLUDE_MULTILEVELANALYZER_HH_

class MultilevelAnalyzer {
public:
	MultilevelAnalyzer(
			int T, int L, int WL_R,
			std::vector<int> level_thickness,
			std::string config_prefix,
			std::vector<int> level_config_num,
			std::vector<void (*)(double*, const double*, int, int, int, int, int, int, int, int)> lowest_level_functions,
			std::vector<std::vector<std::vector<int> > > field_compositions
			);

	std::string config_filename(const std::vector<int>& tag);
	void compute_sublattice_fields(const std::vector<int>& conf_tag, int level, double** T_fields);

private:

	void obtain_sublattice_gauge_field(double*& sub_gauge_field, const std::vector<int>& tag);

	const int T, L, WL_R;
	const std::vector<int> level_thickness;
	const std::string config_prefix;
	const std::vector<int> level_config_num;
	const std::vector<std::vector<std::vector<int> > > field_compositions;
	const std::vector<void (*)(double*, const double*, int, int, int, int, int, int, int, int)> lowest_level_functions;
};

#endif /* INCLUDE_MULTILEVELANALYZER_HH_ */
