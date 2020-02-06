#include <stdexcept>
#include <string>

#include <algorithm/TwolinkComputer.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

TwolinkComputer::TwolinkComputer(std::string name, twolink_operator_sig (*computes_twolink), int t_extent) :
		m_name(name), m_t_extent(t_extent), m_computes_twolink(computes_twolink) {
	if (name.empty() || computes_twolink == nullptr || t_extent < 0)
		throw std::invalid_argument("invalid TwolinkComputer");
}

std::string TwolinkComputer::name() const {
	return m_name;
}

int TwolinkComputer::t_extent() const {
	return m_t_extent;
}

void TwolinkComputer::at(double* result, int t, int x, int y, int z, int dir, int rsep,
		const double* sub_gauge_field, int T, int L) const {
	if (sub_gauge_field == nullptr || T < 1 || L < 1)
		throw std::invalid_argument("invalid gauge field");
	m_computes_twolink(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep);
}

}
}
}
