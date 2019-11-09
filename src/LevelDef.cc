#include <vector>
#include <set>
#include <stdexcept>

#include <TwolinkOperator.hh>
#include <TwolinkOperatorWriter.hh>

#include <LevelDef.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

LevelDef::LevelDef(std::vector<int> timeslice_sizes) :
		m_timeslice_sizes(timeslice_sizes) {

	if (timeslice_sizes.empty())
		throw std::invalid_argument("invalid LevelDef");
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

const std::vector<int>& LevelDef::timeslice_sizes() const {
	return m_timeslice_sizes;
}

const std::vector<TwolinkOperator>& LevelDef::operators() const {
	return m_operators;
}

void LevelDef::add_operator(const TwolinkOperator& def) {
	m_operators.push_back(def);
}

//private

void LevelDef::alloc_operators(const std::set<int>& WL_Rs, int T, int L) {
	for (auto& op : m_operators)
		TwolinkOperatorWriter::alloc_T_fields(op, WL_Rs, m_timeslice_sizes, T, L);
}

void LevelDef::free_operators() {
	for(auto& op : m_operators)
		TwolinkOperatorWriter::free_T_fields(op);
}

}
}
}
