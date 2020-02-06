#include <vector>
#include <map>
#include <set>
#include <utility>
#include <stdexcept>

#include <global_defs.hh>
#include <geometry.hh>

#include <algorithm/T_field.hh>
#include <physics/sublattice_algebra.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

T_field::T_field(const std::vector<std::pair<int, bool> >& thickness_defined, int T, int L) :
		m_T(T), m_L(L) {
	int total_thickness = 0;
	for (const auto& tsl : thickness_defined) {
		if (tsl.first < 1)
			throw std::invalid_argument("timeslice size < 1");
		total_thickness += tsl.first;
	}
	if (T % total_thickness != 0)
		throw std::invalid_argument("timeslices do not fit on the lattice");

	int num_t = 0, curr_t = 0;
	bool first = true;
	while (curr_t < T) {
		for (const auto& tsl : thickness_defined) {
			if (tsl.second) {
				m_t_to_index[curr_t] = num_t++;
				if (first) {
					m_timeslice_size = tsl.first;
					first = false;
				} else if (tsl.first != m_timeslice_size)
					throw std::invalid_argument("cannot define T_field at timeslices of different size");
			}
			curr_t += tsl.first;
		}
	}
	m_data = new double[field_elems()]();
}

T_field::~T_field() {
	delete[] m_data;
}
T_field::T_field(const T_field& fld) {
	m_timeslice_size = fld.m_timeslice_size;
	m_T = fld.m_T;
	m_L = fld.m_L;
	m_t_to_index = fld.m_t_to_index;
	m_data = new double[field_elems()]();
	for (int i = 0; i < field_elems(); ++i)
		m_data[i] = fld.m_data[i];
}
T_field::T_field(T_field&& fld) {
	m_timeslice_size = fld.m_timeslice_size;
	m_T = fld.m_T;
	m_L = fld.m_L;
	m_t_to_index = fld.m_t_to_index;
	m_data = fld.m_data;
	fld.m_data = nullptr;
}
T_field& T_field::operator=(const T_field& fld) {
	m_timeslice_size = fld.m_timeslice_size;
	m_T = fld.m_T;
	m_L = fld.m_L;
	m_t_to_index = fld.m_t_to_index;
	m_data = new double[field_elems()]();
	for (int i = 0; i < field_elems(); ++i)
		m_data[i] = fld.m_data[i];
	return *this;
}
T_field& T_field::operator=(T_field&& fld) {
	m_timeslice_size = fld.m_timeslice_size;
	m_T = fld.m_T;
	m_L = fld.m_L;
	m_t_to_index = fld.m_t_to_index;
	m_data = fld.m_data;
	fld.m_data = nullptr;
	return *this;
}

std::set<int> T_field::defined_ts() const {
	std::set<int> ts;
	for (auto t_i : m_t_to_index)
		ts.insert(t_i.first);
	return ts;
}

double* T_field::T_at(int t, int x, int y, int z, int dir) const {
	return m_data + (3 * get_index(m_t_to_index.at((t + m_T) % m_T), x, y, z, m_t_to_index.size(), m_L) + (dir - 1)) * SO_elems;
}

T_field& T_field::operator/=(double d) {
	for (int i = 0; i < field_elems(); ++i)
		m_data[i] /= d;
	return *this;
}

}
}
}
