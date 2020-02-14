#include <iostream>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <numeric>

#include <helper_functions.hh>
#include <io_tools.hh>

namespace de_uni_frankfurt_itp::reisinger::multilevel_0819::find_sizes {

using t_coord = int;
using t_extent = int;
using timeslice_sizes = std::vector<int>;

enum observable_mode {
	center
};

void print_syntax_help(char* argv0) {
	std::cout << "\nUsage: " << argv0 << ""
			"\t[--help] [--require <t_extents>] [--sizes <pattern_sizes>] <tsl_num> <min_size> <max_size> <t_extents>\n\n";
}

void print_option_help() {
	std::cout
			<< "\n"
					"Search all timeslice patterns consisting of <tsl_num> timeslices with sizes in the range [<min_size>,<max_size>]\n"
					"for possible placements of operators with temporal extent in the list <t_extents>. Timeslice patterns are repeated\n"
					"periodically as needed. Operators can be placed if no timeslice boundary lies at a set of time coordinates\n"
					"determined by <mode>. For each pattern, prints possible placements:\n"
					"\t[<total_pattern_size>] { <pattern> } : <extent_1> { <time_coord_1_1> <time_coord_1_2> ... } <extent_2> { <time_coord_2_1> ... } ...\n"
					"\n"
					"Parameters\n"
					"\n"
					"\t<mode>\n"
					"\t\tNot implemented (always default = center). Choose a set of time coordinates relative to the operator position\n"
					"\t\tthat cannot lie on a timeslice boundary.\n"
					"\n"
					"\n"
					"Options\n"
					"\n"
					"\t--require | -r <t_extents>\n"
					"\t\tOnly print result for a given pattern if at least all of the extents in the list <t_extents> are possible\n"
					"\n"
					"\t--sizes | -s <pattern_sizes>\n"
					"\t\tOnly search patterns with total size in the list <pattern_sizes>"
					"\n"
					"\n"
					"Modes\n"
					"\n"
					"\tcenter\n"
					"\t\tThe midpoint of the temporal Wilson line. For odd temporal extent, shifted by +/- half a lattice site.\n"
					"\n";
}

bool handle_GNU_options(int argc, char**& argv, std::set<t_extent>& required_extents, std::set<int>& total_pattern_sizes) {
	constexpr int REQUIRED_ARG_NUM = 4;

	static struct option long_opts[] = {
			{ "help", no_argument, 0, 'h' },
			{ "require", required_argument, 0, 'r' },
			{ "sizes", required_argument, 0, 's' },
			{ 0, 0, 0, 0 }
	};

	int opt = -1, long_opts_i = 0;
	while ((opt = getopt_long(argc, argv, "hr:s:", long_opts, &long_opts_i)) != -1) {
		switch (opt) {
			case 'h':
				print_syntax_help(argv[0]);
				print_option_help();
				return false;
				break;
			case 'r': {
				std::vector<t_extent> required_extents_vec = tools::helper::parse_unsigned_int_list(optarg);
				required_extents = std::set<t_extent>(required_extents_vec.begin(), required_extents_vec.end());
				break;
			}
			case 's': {
				std::vector<int> pattern_sizes_vec = tools::helper::parse_unsigned_int_list(optarg);
				total_pattern_sizes = std::set<int>(pattern_sizes_vec.begin(), pattern_sizes_vec.end());
				break;
			}
		}
	}

	if (argc - optind != REQUIRED_ARG_NUM) {
		print_syntax_help(argv[0]);
		return false;
	}

	argv = argv + optind - 1;

	return true;
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
		const timeslice_sizes& timeslice_sizes_pattern, const std::set<std::set<t_coord> >& boundary_forbidden_time_coords) {

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
		case center:
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
		const std::set<t_extent>& t_extents, observable_mode mode, const std::set<int>& total_pattern_sizes,
		std::map<timeslice_sizes, std::map<t_extent, std::set<t_coord> > >& possible_placements) {

	if (total_pattern_sizes.empty() || total_pattern_sizes.count(
			std::accumulate(timeslice_sizes_pattern.begin(), timeslice_sizes_pattern.end(), 0))) {
		auto possible_extents = try_place_operator(timeslice_sizes_pattern, t_extents, mode);
		if (!possible_extents.empty())
			possible_placements[timeslice_sizes_pattern] = possible_extents;
	}
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

timeslice_sizes cycle_sizes(const timeslice_sizes& sizes) {

	timeslice_sizes cycled_sizes(sizes);

	if (sizes.size() < 2)
		return sizes;

	cycled_sizes.insert(cycled_sizes.begin(), cycled_sizes.back());
	cycled_sizes.pop_back();

	return cycled_sizes;
}

void remove_cyclic(std::map<timeslice_sizes, std::map<t_extent, std::set<t_coord> > >& operator_placements) {
	for (const auto& [sizes, extents] : operator_placements) {
		(void) extents;
		for (timeslice_sizes equal_sizes = cycle_sizes(sizes); equal_sizes != sizes; equal_sizes = cycle_sizes(equal_sizes))
			if (auto eq_pos = operator_placements.find(equal_sizes); eq_pos != operator_placements.end())
				operator_placements.erase(eq_pos);
	}
}

void remove_reverse(std::map<timeslice_sizes, std::map<t_extent, std::set<t_coord> > >& operator_placements) {
	for (const auto& [sizes, extents] : operator_placements) {
		(void) extents;
		auto reverse = sizes;
		std::reverse(reverse.begin(), reverse.end());
		if (reverse != sizes)
			if (auto rev_pos = operator_placements.find(reverse); rev_pos != operator_placements.end())
				operator_placements.erase(rev_pos);
	}
}

void remove_equal_timeslice_patterns(std::map<timeslice_sizes, std::map<t_extent, std::set<t_coord> > >& operator_placements) {
	remove_reverse(operator_placements);
	remove_cyclic(operator_placements);
}

void print(const std::map<timeslice_sizes, std::map<t_extent, std::set<t_coord> > >& operator_placements) {
	using tools::io_tools::operator<<;

	for (const auto& pattern : operator_placements) {
		std::cout << "[" << std::accumulate(pattern.first.begin(), pattern.first.end(), 0) << "] " << pattern.first << " : ";
		for (const auto& t_extent : pattern.second)
			std::cout << t_extent.first << " " << t_extent.second << " ";
		std::cout << "\n";
	}
}

}

int main(int argc, char** argv) {
	using namespace std;
	using namespace de_uni_frankfurt_itp::reisinger;
	using namespace multilevel_0819::find_sizes;
	using tools::helper::parse_unsigned_int_list;
	using tools::helper::nest_for;

	set<t_extent> required_extents;
	set<int> total_pattern_sizes;
	if (!handle_GNU_options(argc, argv, required_extents, total_pattern_sizes))
		return 0;

	int arg_num = 1;
	const int tsl_num = stoi(argv[arg_num++]);
	const int min_size = stoi(argv[arg_num++]);
	const int max_size = stoi(argv[arg_num++]);

	const vector<t_extent> t_extents_vec = parse_unsigned_int_list(argv[arg_num]);
	const set<t_extent> t_extents = set<t_extent>(t_extents_vec.begin(), t_extents_vec.end());

//	*********************************************************************************************************

	const vector<pair<int, int> > size_limits(tsl_num, { min_size, max_size + 1 });

	map<timeslice_sizes, map<t_extent, set<t_coord> > > possible_operator_placements;
	nest_for(size_limits, scan_operator_placements,
			t_extents, center, total_pattern_sizes, possible_operator_placements);

	remove_equal_timeslice_patterns(possible_operator_placements);
	filter(possible_operator_placements, required_extents);

	print(possible_operator_placements);

}
