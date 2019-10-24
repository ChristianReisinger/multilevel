#include <string>

#include <twolink_operator_functions.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_OPERATORFACTOR_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_OPERATORFACTOR_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class OperatorFactor {
public:
	virtual std::string name() const = 0;
	virtual int t_extent() const = 0;
	virtual void at(double* result, int t, int x, int y, int z, int dir, int rsep,
			const double* sub_gauge_field = nullptr, int T = 0, int L = 0) const = 0;

protected:
	virtual ~OperatorFactor() = default;
};

}
}
}

#endif
