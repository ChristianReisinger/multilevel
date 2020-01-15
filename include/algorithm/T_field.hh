#include <vector>
#include <set>
#include <utility>
#include <map>

#include <physics/sublattice_algebra.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_T_FIELD_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_T_FIELD_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class T_field {
public:
	T_field(const std::vector<std::pair<int, bool> >& thickness_defined, int T, int L);

	~T_field();
	T_field(const T_field& fld);
	T_field(T_field&& fld);
	T_field& operator=(const T_field& fld);
	T_field& operator=(T_field&& fld);

	std::set<int> defined_ts() const;
	double* T_at(int t, int x, int y, int z, int dir) const;
	inline int timeslice_size() const {
		return m_timeslice_size;
	}
	T_field& operator/=(double d);

private:

	inline int field_elems() {
		return m_t_to_index.size() * m_L * m_L * m_L * 3 * SO_elems;
	}

	double* m_data;
	int m_timeslice_size, m_T, m_L;
	std::map<int, int> m_t_to_index;
};

}
}
}

#endif
