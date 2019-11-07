#include <vector>
#include <stdexcept>
#include <utility>
#include <algorithm>

#include <global_defs.hh>
#include <linear_algebra.hh>

#include <FactorInterface.hh>
#include <T_field.hh>
#include <sublattice_algebra.hh>

#include <TwolinkOperator.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

TwolinkOperator::TwolinkOperator(std::string name, std::vector<bool> timeslice_isdefined,
		std::vector<const FactorInterface*> factors) :
		m_name(name), m_timeslice_isdefined(timeslice_isdefined), factors(factors), m_t_extent(0) {

	if (name.empty() || factors.empty() || !valid_timeslice_def())
		throw std::invalid_argument("invalid TwolinkOperator");

	for (const FactorInterface* def : factors)
		m_t_extent += def->t_extent();

}

int TwolinkOperator::timeslice_num_per_cycle() const {
	return std::count_if(m_timeslice_isdefined.begin(), m_timeslice_isdefined.end(), [](bool b) {return b;});
}

std::string TwolinkOperator::descr() const {
	return m_descr;
}

std::string TwolinkOperator::descr(std::string str) {
	m_descr = str;
	return descr();
}
std::string TwolinkOperator::name() const {
	return m_name;
}

int TwolinkOperator::t_extent() const {
	return m_t_extent;
}

void TwolinkOperator::at(double* result, int t, int x, int y, int z, int dir, int rsep,
		const double* sub_gauge_field, int T, int L) const {

	so_eq_so(result, m_r_fields.at(rsep).T_at(t, x, y, z, dir));
}

int TwolinkOperator::timeslice_num() const {
	int tsl_num = 0;
	for (bool b : m_timeslice_isdefined)
		if (b)
			++tsl_num;
	return tsl_num;
}

std::set<int> TwolinkOperator::defined_ts(int WL_R) const {
	return m_r_fields.at(WL_R).defined_ts();
}

//private

void TwolinkOperator::alloc_T_fields(const std::set<int>& WL_Rs,
		const std::vector<int>& timeslice_sizes, int T, int L) {

	std::vector<std::pair<int, bool> > timeslice_size_defined;
	if (timeslice_sizes.size() != m_timeslice_isdefined.size())
		throw std::invalid_argument("invalid timeslice definitions");
	for (int tsl = 0; tsl < timeslice_sizes.size(); ++tsl)
		timeslice_size_defined.emplace_back(timeslice_sizes[tsl], m_timeslice_isdefined[tsl]);

	for (int WL_R : WL_Rs)
		m_r_fields.emplace(WL_R, T_field(timeslice_size_defined, T, L));
}

void TwolinkOperator::free_T_fields() {
	m_r_fields.clear();
}

bool TwolinkOperator::valid_timeslice_def() {
	bool valid = false;
	for (bool b : m_timeslice_isdefined)
		if (b) {
			valid = true;
			break;
		}
	return valid;
}

}
}
}
