/*
 * sublattice_fields.hh
 *
 *  Created on: 9 Aug 2019
 *      Author: reisinger
 */

#include <vector>
#include <set>
#include <utility>
#include <map>

#include <geometry.hh>
#include <sublattice_algebra.hh>

#ifndef INCLUDE_T_FIELD_HH_
#define INCLUDE_T_FIELD_HH_

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
	inline int timeslice_thickness() const {
		return thickness;
	}
	T_field& operator/=(double d);

private:

	inline int field_elems() {
		return t_to_index.size() * L * L * L * 3 * SO_elems;
	}

	double* data;
	int thickness, T, L;
	std::map<int, int> t_to_index;
};

}
}
}

#endif /* INCLUDE_T_FIELD_HH_ */
