#include <set>
#include <vector>

#include <TwolinkOperator.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_TWOLINKOPERATORWRITER_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_TWOLINKOPERATORWRITER_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class TwolinkOperatorWriter {
private:
	inline static void alloc_T_fields(TwolinkOperator& op, const std::set<int>& WL_Rs,
			const std::vector<int>& timeslice_sizes, int T, int L) {
		op.alloc_T_fields(WL_Rs, timeslice_sizes, T, L);
	}
	inline static T_field& field(TwolinkOperator& op, int WL_R) {
		return op.m_r_fields.at(WL_R);
	}
	friend class LevelDef;
	friend class MultilevelAnalyzer;
};

}
}
}

#endif
