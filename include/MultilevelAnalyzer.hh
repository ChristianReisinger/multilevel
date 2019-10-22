#include <vector>
#include <set>
#include <map>
#include <string>
#include <chrono>

#include <MultilevelConfig.hh>

#ifndef INCLUDE_MULTILEVELANALYZER_HH_
#define INCLUDE_MULTILEVELANALYZER_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class MultilevelAnalyzer {
public:
	MultilevelAnalyzer(MultilevelConfig& multilevel_config, std::set<int> WL_Rs,
			std::vector<std::map<std::string, std::vector<std::string> > > level_operator_factors,
			std::vector<std::map<std::string, std::vector<bool> > > level_operator_timeslice_defined,
			std::map<std::string, void (*)(double*, const double*, int, int, int&, int, int, int, int, int)> lowest_level_functions
			);

	std::map<std::string, std::map<int, T_field> > compute_T_fields();
	int milliseconds_spent_computing();

private:

	void alloc_T_fields(std::map<std::string, std::map<int, T_field> >& T_fields, const int level);
	void compute_sublattice_fields(std::map<std::string, std::map<int, T_field> >& T_fields, const int level);

	MultilevelConfig* config;
	const std::set<int> WL_Rs;
	const std::vector<std::map<std::string, std::vector<std::string> > > level_operator_factors;
	const std::vector<std::map<std::string, std::vector<bool> > > level_operator_timeslice_defined;
	const std::map<std::string, void (*)(double*, const double*, int, int, int&, int, int, int, int, int)> lowest_level_functions;

	std::chrono::steady_clock::duration time_spent_computing_operators;
};

}
}
}

#endif /* INCLUDE_MULTILEVELANALYZER_HH_ */
