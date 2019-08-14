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
	/**
	 * @param generate_configs	if false, read existing configs with tags {i,j,..,k} for levels 0,1,.. and filenames ending in .i.j...k (with zero padding)
	 * 							otherwise, generate configs successively during the algorithm, starting from a thermalized one with tag { 1 }
	 *
	 *
	 * @param field_compositions a vector V at index L in field_compositions describes compositions
	 *							of the fields F_i(L) at level L from those at level L+1.
	 *							A vector W at index i in V contains indices j of the fields at level L+1.
	 *							The field F_i(L) is obtained by multiplying fields F_j(L+1) for all j in the order
	 *							they appear in W.
	 *							At the lowest level (= field_compositions.size() - 1), j instead is an index of a
	 *							function in lowest_level_functions.
	 */
	MultilevelAnalyzer(
			int T, int L, int WL_R,
			std::vector<int> level_thickness,
			std::string config_prefix,
			std::vector<int> level_config_num,
			std::vector<void (*)(double*, const double*, int, int, int, int, int, int, int, int)> lowest_level_functions,
			std::vector<std::vector<std::vector<int> > > field_compositions,
			bool generate_configs = false, double beta = 0.0, int seed = 0, std::vector<int> level_updates = {}
			);
	~MultilevelAnalyzer();

	std::string tag_to_string(const std::vector<int>& tag);
	std::string config_filename(const std::vector<int>& tag);
	void compute_sublattice_fields(const std::vector<int>& conf_tag, int level, double** T_fields);

private:

	void read_sublattice_gauge_field(double*& sub_gauge_field, const std::vector<int>& tag);
	void update_sublattice_gauge_field(const std::vector<int>& tag);

	const bool generate_configs;
	const double beta;
	const int seed;
	const std::vector<int> level_updates;
	double* config_buf;

	const int T, L, WL_R;
	const std::vector<int> level_thickness;
	const std::string config_prefix;
	const std::vector<int> level_config_num;
	const std::vector<std::vector<std::vector<int> > > field_compositions;
	const std::vector<void (*)(double*, const double*, int, int, int, int, int, int, int, int)> lowest_level_functions;
};

#endif /* INCLUDE_MULTILEVELANALYZER_HH_ */
