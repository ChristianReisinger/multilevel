#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>

#include <OperatorFactor.hh>
#include <T_field.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_TWOLINKOPERATOR_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_TWOLINKOPERATOR_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class TwolinkOperator: OperatorFactor {
public:
	TwolinkOperator(std::string name, std::vector<bool> timeslice_isdefined, std::vector<OperatorFactor*> factors);

	int timeslice_num_per_cycle() const;
	std::string descr() const;
	std::string descr(std::string str);

	std::string name() const override;
	int t_extent() const override;
	void at(double* result, int t, int x, int y, int z, int dir, int rsep,
			const double* sub_gauge_field = nullptr, int T = 0, int L = 0) const override;

	void alloc_T_fields(const std::set<int>& WL_Rs, const std::vector<int>& timeslice_sizes, int T, int L);

private:
	std::string m_name, m_descr = "";
	std::vector<bool> m_timeslice_isdefined;
	std::vector<OperatorFactor*> m_factors;

	std::string m_name;
	int m_t_extent;
	std::map<int, T_field> m_r_fields { };

	bool valid_timeslice_def();

};

}
}
}

#endif
