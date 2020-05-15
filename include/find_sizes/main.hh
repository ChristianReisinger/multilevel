#include <vector>
#include <utility>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_FIND_SIZES_MAIN_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_FIND_SIZES_MAIN_HH_

namespace de_uni_frankfurt_itp::reisinger::multilevel_0819::find_sizes {

using t_coord = int;
using t_extent = int;
using timeslice_sizes = std::vector<int>;
using timeslice_setups = std::map<timeslice_sizes, std::map<t_extent, std::set<t_coord> > >;

enum observable_mode {
	center
};

timeslice_setups find_sizes(const std::vector<std::pair<int, int> >& size_limits, const std::set<int>& total_pattern_sizes,
		const std::set<t_extent>& t_extents, observable_mode mode);

void filter(timeslice_setups& setups, const std::set<t_extent>& required_extents);

void print(const timeslice_setups& setups);

}

#endif
