#include <vector>
#include <set>
#include <map>
#include <string>
#include <chrono>

#include <LevelDef.hh>

#include <MultilevelConfig.hh>

#ifndef INCLUDE_MULTILEVELANALYZER_HH_
#define INCLUDE_MULTILEVELANALYZER_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class MultilevelAnalyzer {
public:
	MultilevelAnalyzer(std::vector<LevelDef>& levels, MultilevelConfig& multilevel_config, std::set<int> WL_Rs);

	std::map<std::string, std::map<int, T_field> > compute_T_fields();
	int milliseconds_spent_computing();

private:

	void alloc_T_fields(std::map<std::string, std::map<int, T_field> >& T_fields, const int level);
	void compute_sublattice_fields(std::map<std::string, std::map<int, T_field> >& T_fields, const int level);

	MultilevelConfig* config;
	const std::set<int> WL_Rs;
	std::vector<LevelDef*> m_levels;

	std::chrono::steady_clock::duration time_spent_computing_operators { 0 };
};

}
}
}

#endif /* INCLUDE_MULTILEVELANALYZER_HH_ */
