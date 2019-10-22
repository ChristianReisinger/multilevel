#include <vector>
#include <map>
#include <set>
#include <utility>
#include <stdexcept>

#include <global_defs.hh>
#include <geometry.hh>
#include <sublattice_algebra.hh>

#include <T_field.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

T_field::T_field(const std::vector<std::pair<int, bool> >& thickness_defined, int T, int L) :
		T(T), L(L) {
	int total_thickness = 0;
	for (const auto& tsl : thickness_defined) {
		if (tsl.first < 1)
			throw std::invalid_argument("timeslice thickness < 1");
		total_thickness += tsl.first;
	}
	if (T % total_thickness != 0)
		throw std::invalid_argument("timeslices do not fit on the lattice");

	int num_t = 0, curr_t = 0;
	bool first = true;
	while (curr_t < T) {
		for (const auto& tsl : thickness_defined) {
			if (tsl.second) {
				t_to_index[curr_t] = num_t++;
				if (first) {
					thickness = tsl.first;
					first = false;
				} else if (tsl.first != thickness)
					throw std::invalid_argument("cannot define T_field at different thickness");
			}
			curr_t += tsl.first;
		}
	}
	data = new double[field_elems()]();
}

T_field::~T_field() {
	delete[] data;
}
T_field::T_field(const T_field& fld) {
	thickness = fld.thickness;
	T = fld.T;
	L = fld.L;
	t_to_index = fld.t_to_index;
	data = new double[field_elems()]();
	for (int i = 0; i < field_elems(); ++i)
		data[i] = fld.data[i];
}
T_field::T_field(T_field&& fld) {
	thickness = fld.thickness;
	T = fld.T;
	L = fld.L;
	t_to_index = fld.t_to_index;
	data = fld.data;
	fld.data = nullptr;
}
T_field& T_field::operator=(const T_field& fld) {
	thickness = fld.thickness;
	T = fld.T;
	L = fld.L;
	t_to_index = fld.t_to_index;
	data = new double[field_elems()]();
	for (int i = 0; i < field_elems(); ++i)
		data[i] = fld.data[i];
	return *this;
}
T_field& T_field::operator=(T_field&& fld) {
	thickness = fld.thickness;
	T = fld.T;
	L = fld.L;
	t_to_index = fld.t_to_index;
	data = fld.data;
	fld.data = nullptr;
	return *this;
}

std::set<int> T_field::defined_ts() const {
	std::set<int> ts;
	for (auto t_i : t_to_index)
		ts.insert(t_i.first);
	return ts;
}

double* T_field::T_at(int t, int x, int y, int z, int dir) const {
	return data + (3 * get_index(t_to_index.at((t + T) % T), x, y, z, t_to_index.size(), L) + (dir - 1)) * SO_elems;
}

T_field& T_field::operator/=(double d) {
	for (int i = 0; i < field_elems(); ++i)
		data[i] /= d;
	return *this;
}

}
}
}
