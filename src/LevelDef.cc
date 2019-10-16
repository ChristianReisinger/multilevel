#include <vector>
#include <stdexcept>

#include <LevelDef.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

LevelDef::LevelDef(std::vector<int> timeslice_sizes) :
		m_timeslice_sizes(timeslice_sizes) {

	if (timeslice_sizes.empty())
		throw std::invalid_argument("invalid level");
	for (const int tsl : timeslice_sizes)
		if (tsl < 1)
			throw std::invalid_argument("invalid timeslice size");
}

int LevelDef::config_num() const {
	return m_config_num;
}
int LevelDef::config_num(int configs) {
	if (configs < 1)
		throw std::invalid_argument("invalid number of configs");
	m_config_num = configs;
	return config_num();
}
int LevelDef::update_num() const {
	return m_update_num;
}
int LevelDef::update_num(int updates) {
	if (updates < 0)
		throw std::invalid_argument("invalid number of updates");
	m_update_num = updates;
	return update_num();
}

const std::vector<TwolinkOperator>& LevelDef::operators() const {
	return m_operators;
}

void LevelDef::add_operator(const TwolinkOperator& def) {
	//TODO check for collision
	m_operators.push_back(def);
}

void LevelDef::alloc_operators(const std::set<int>& WL_Rs, int T, int L) {
	for (auto& op : m_operators) {
		op.alloc_T_fields(WL_Rs, m_timeslice_sizes, T, L);
	}
}

}
}
}
