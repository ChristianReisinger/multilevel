#include <string>
#include <stdexcept>

#include <OperatorFactor.hh>
#include <twolink_operator_functions.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_TWOLINKCOMPUTER_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_TWOLINKCOMPUTER_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class TwolinkComputer: public OperatorFactor {
public:
	TwolinkComputer(std::string name, twolink_operator_sig (*computes_twolink), int t_extent);

	std::string name() const override;
	int t_extent() const override;
	void at(double* result, int t, int x, int y, int z, int dir, int rsep,
			const double* sub_gauge_field, int T, int L) const override;

private:
	std::string m_name;
	int m_t_extent;
	twolink_operator_sig (*m_computes_twolink);
};

}
}
}

#endif
