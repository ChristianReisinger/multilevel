/*
 * MultilevelAnalyzer.hh
 *
 *  Created on: 9 Aug 2019
 *      Author: reisinger
 */

#include <vector>
#include <string>
#include <chrono>

#include <MultilevelConfig.hh>

#ifndef INCLUDE_MULTILEVELANALYZER_HH_
#define INCLUDE_MULTILEVELANALYZER_HH_

class MultilevelAnalyzer {
public:
	/**
	 * @param field_compositions a vector V at index L in field_compositions describes compositions
	 *							of the fields F_i(L) at level L from those at level L+1.
	 *							A vector W at index i in V contains indices j of the fields at level L+1.
	 *							The field F_i(L) is obtained by multiplying fields F_j(L+1) for all j in the order
	 *							they appear in W.
	 *							At the lowest level (= field_compositions.size() - 1), j instead is an index of a
	 *							function in lowest_level_functions.
	 */
	MultilevelAnalyzer(MultilevelConfig& multilevel_config, int WL_R,
			std::vector<std::vector<std::vector<int> > > field_compositions,
			std::vector<void (*)(double*, const double*, int, int, int&, int, int, int, int, int)> lowest_level_functions
			);

	void compute_sublattice_fields(double** T_fields, const int level = 0);
	int milliseconds_spent_computing();

private:

	MultilevelConfig* config;
	const int WL_R;
	const std::vector<std::vector<std::vector<int> > > field_compositions;
	const std::vector<void (*)(double*, const double*, int, int, int&, int, int, int, int, int)> lowest_level_functions;

	std::chrono::steady_clock::duration time_spent_computing_operators;

	bool top_level = true;
};

#endif /* INCLUDE_MULTILEVELANALYZER_HH_ */
