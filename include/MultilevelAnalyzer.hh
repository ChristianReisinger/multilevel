#include <vector>
#include <set>
#include <map>
#include <string>
#include <chrono>

#include <LevelDef.hh>

#include <MultilevelConfig.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_MULTILEVELANALYZER_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_MULTILEVELANALYZER_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class MultilevelAnalyzer {
public:
	MultilevelAnalyzer(std::vector<LevelDef>& levels, MultilevelConfig& multilevel_config, std::set<int> WL_Rs);
	~MultilevelAnalyzer() = default;
	MultilevelAnalyzer(const MultilevelAnalyzer&) = delete;
	MultilevelAnalyzer(MultilevelAnalyzer&&) = delete;
	MultilevelAnalyzer& operator=(const MultilevelAnalyzer&) = delete;
	MultilevelAnalyzer& operator=(MultilevelAnalyzer&&) = delete;

	void compute_T_fields();
	int milliseconds_spent_computing() const;

private:

	void compute_sublattice_fields(const size_t level);

	std::vector<LevelDef*> m_levels;
	MultilevelConfig* const m_config;
	const std::set<int> m_WL_Rs;

	std::chrono::steady_clock::duration time_spent_computing_operators { 0 };
};

}
}
}

#endif
