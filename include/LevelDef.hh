#include <vector>
#include <set>

#include <TwolinkOperator.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_LEVELDEF_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_LEVELDEF_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class LevelDef {
public:
	LevelDef(std::vector<int> timeslice_sizes, int T, int L);

	int config_num() const;
	int config_num(int configs);
	int update_num() const;
	int update_num(int updates);
	const std::vector<int>& timeslice_sizes() const;

	const std::vector<TwolinkOperator>& operators() const;
	void add_operator(const TwolinkOperator& def);
	void alloc_operators(const std::set<int>& WL_Rs);

private:
	const int m_T, m_L;
	int m_config_num = 1, m_update_num = 0;
	std::vector<int> m_timeslice_sizes;
	std::vector<TwolinkOperator> m_operators { };
};

}
}
}

#endif
