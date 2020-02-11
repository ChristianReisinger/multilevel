#include <iostream>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <numeric>

#include <helper_functions.hh>
#include <io_tools.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {
namespace try_operator_on_timeslices {

using t_coord = int;
using t_extent = int;
using timeslice_sizes = std::vector<int>;

enum observable_mode {
	force
};

void print_syntax_help(char* argv0) {
	std::cout << "\nUsage: " << argv0 << ""
			"\t[--help] [-b <no_boundary_times>] <timeslice_sizes> <t_extent>\n\n";
}

void print_option_help() {
	std::cout << "\n"
			"Parameters\n"
			"\n"
			"\t<no_boundary_times>\n"
			"\t\tlist of times relative to the temporal coordinate of the operator\n"
			"\t\twhich cannot lie on a timeslice boundary\n"
			"\n"
			"\t<timeslice_sizes>\n"
			"\t\tlist of timeslize sizes. If the total size is smaller than the operator size,\n"
			"\t\tthe pattern is repeated periodically.\n"
			"\n"
			"\t<t_extent>\n"
			"\t\ttemporal operator extent\n"
			"\n";
}

bool handle_GNU_options(int argc, char**& argv, std::vector<int>& no_boundary_times) {
	const int REQUIRED_ARG_NUM = 2;
	bool no_help_required = true;

	static struct option long_opts[] = {
			{ "help", no_argument, 0, 'h' },
			{ "noboundary", required_argument, 0, 'b' },
			{ 0, 0, 0, 0 }
	};

	int opt = -1, long_opts_i = 0;
	while ((opt = getopt_long(argc, argv, "hb:", long_opts, &long_opts_i)) != -1) {
		switch (opt) {
			case 'h':
				print_syntax_help(argv[0]);
				print_option_help();
				no_help_required = false;
				break;
			case 'b':
				no_boundary_times = de_uni_frankfurt_itp::reisinger::tools::helper::parse_unsigned_int_list(optarg);
				break;
		}
	}

	if (no_help_required && argc - optind != REQUIRED_ARG_NUM) {
		print_syntax_help(argv[0]);
		no_help_required = false;
	}

	argv = argv + optind - 1;

	return no_help_required;
}

int next_timeslice_boundary(int t, const timeslice_sizes& timeslice_sizes_pattern) {
	const int tsl_num = timeslice_sizes_pattern.size();
	int curr_t = 0;

	for (int tsl_i = 0; curr_t <= t; ++tsl_i)
		curr_t += timeslice_sizes_pattern[tsl_i % tsl_num];
	return curr_t;
}

bool is_timeslice_boundary(int t, const timeslice_sizes& timeslice_sizes_pattern) {
	int curr_t = 0;
	while (curr_t <= t) {
		if (t == curr_t)
			return true;

		curr_t = next_timeslice_boundary(curr_t, timeslice_sizes_pattern);
	}
	return false;
}

/**
 * @param operator_t_extent					size of the operator in temporal direction
 * @param boundary_forbidden_time_coords	for each set in no_boundary_times, at least of the times in that set must
 * 											not be a timeslice boundary (i.e. elements in the outer set are related by
 * 											&&, elements in the inner sets by ||)
 * @return	set of time coordinates where an operator can be placed to start and end on a timeslice
 * 			boundary, and fulfill the required conditions given by the parameters
 */
auto find_possible_operator_time_coords(t_extent operator_t_extent,
		const timeslice_sizes& timeslice_sizes_pattern, const std::set<std::set<int> >& boundary_forbidden_time_coords) {

	std::set<t_coord> possible_times;

	for (int operator_pos_t = 0; operator_pos_t < accumulate(timeslice_sizes_pattern.begin(), timeslice_sizes_pattern.end(), 0);
			operator_pos_t = next_timeslice_boundary(operator_pos_t, timeslice_sizes_pattern)) {

		bool non_boundaries_fulfilled = true;
		for (const auto& one_of_t_not_a_boundary : boundary_forbidden_time_coords) {
			bool found_non_boundary = false;
			for (const int t : one_of_t_not_a_boundary) {
				if (!is_timeslice_boundary(operator_pos_t + t, timeslice_sizes_pattern)) {
					found_non_boundary = true;
					break;
				}
			}
			if (!found_non_boundary) {
				non_boundaries_fulfilled = false;
				break;
			}
		}

		if (non_boundaries_fulfilled && is_timeslice_boundary(operator_pos_t + operator_t_extent, timeslice_sizes_pattern))
			possible_times.insert(operator_pos_t);
	}

	return possible_times;
}

auto get_boundary_forbidden_time_coords(t_extent operator_t_extent, observable_mode mode) {

	std::set<std::set<t_coord> > boundary_forbidden_t_coords;

	switch (mode) {
		case force:
			std::set<t_extent> center_coords = { operator_t_extent / 2 };
			if (operator_t_extent % 2 != 0)
				center_coords.insert(operator_t_extent / 2 + 1);
			boundary_forbidden_t_coords.insert(center_coords);
			break;
	}
	return boundary_forbidden_t_coords;
}

auto try_place_operator(const timeslice_sizes& timeslice_sizes_pattern,
		const std::set<t_extent>& t_extents, observable_mode mode) {

	std::map<t_extent, std::set<t_coord> > possible_extents;

	for (const t_extent operator_t_extent : t_extents) {
		std::set<t_coord> possible_operator_t_coords = find_possible_operator_time_coords(operator_t_extent,
				timeslice_sizes_pattern, get_boundary_forbidden_time_coords(operator_t_extent, mode));

		if (!possible_operator_t_coords.empty())
			possible_extents[operator_t_extent] = possible_operator_t_coords;
	}

	return possible_extents;
}

void scan_operator_placements(const timeslice_sizes& timeslice_sizes_pattern,
		const std::set<t_extent>& t_extents, observable_mode mode,
		std::map<timeslice_sizes, std::map<t_extent, std::set<t_coord> > >& possible_placements) {

	auto possible_extents = try_place_operator(timeslice_sizes_pattern, t_extents, mode);
	if (!possible_extents.empty())
		possible_placements[timeslice_sizes_pattern] = possible_extents;
}

void filter(std::map<timeslice_sizes, std::map<t_extent, std::set<t_coord> > >& operator_placements,
		const std::set<t_extent>& required_extents) {

	for (auto it = operator_placements.begin(); it != operator_placements.end();) {
		bool has_required = true;
		for (t_extent ext : required_extents)
			if (!it->second.count(ext)) {
				has_required = false;
				break;
			}
		if (has_required)
			++it;
		else
			it = operator_placements.erase(it);
	}
}

void print(const std::map<timeslice_sizes, std::map<t_extent, std::set<t_coord> > >& operator_placements) {
	using tools::io_tools::operator<<;

	for (const auto& pattern : operator_placements) {
		std::cout << pattern.first << " : ";
		for (const auto& t_extent : pattern.second)
			std::cout << t_extent.first << " " << t_extent.second << " ";
		std::cout << "\n";
	}
}

}
}
}
}

int main(int argc, char** argv) {
	using namespace std;
	using namespace de_uni_frankfurt_itp::reisinger;
	using namespace multilevel_0819::try_operator_on_timeslices;
	using tools::helper::parse_unsigned_int_list;
	using tools::helper::nest_for;

	/*TODO change program syntax, parse new parameters
	 * options: -f <required_extents>
	 */
	const int tsl_num = 4, min_size = 1, max_size = 3;
	const set<t_extent> t_extents = { 4, 5, 6 };
	const vector<pair<int, int> > size_limits(tsl_num, { min_size, max_size + 1 });

	map<timeslice_sizes, map<t_extent, set<t_coord> > > possible_operator_placements;
	nest_for(size_limits,
			scan_operator_placements, t_extents, force, possible_operator_placements);

	filter(possible_operator_placements, t_extents);

	print(possible_operator_placements);

}
