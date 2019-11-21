#include <set>
#include <vector>

#include <LevelDef.hh>
#include <TwolinkOperator.hh>
#include <MultilevelConfig.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_LEVELACCESS_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_LEVELACCESS_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class LevelAccess {
private:
	inline static std::vector<TwolinkOperator>& operators(LevelDef& level) {
		return level.m_operators;
	}
	inline static void alloc_operators(LevelDef& level, const std::set<int>& WL_Rs, int T, int L) {
		level.alloc_operators(WL_Rs, T, L);
	}
	inline static void free_operators(LevelDef& level) {
		level.free_operators();
	}
	inline static void set_levels(MultilevelConfig& config, std::vector<LevelDef*> levels) {
		config.set_levels(levels);
	}
	inline static void update(MultilevelConfig& config, size_t level) {
		config.update(level);
	}

	friend class MultilevelAnalyzer;
};

}
}
}

#endif
