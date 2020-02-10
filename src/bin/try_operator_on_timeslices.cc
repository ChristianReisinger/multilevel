#include <iostream>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <numeric>

#include <helper_functions.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {
namespace try_operator_on_timeslices {

using t_coord = int;
using t_extent = int;
using timeslice_sizes = std::vector<int>;

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
 * @param operator_t_extent	size of the operator in temporal direction
 * @param no_boundary_times	for each set in no_boundary_times, at least of the times in that set must
 * 							not be a timeslice boundary (i.e. elements in the outer set are related by
 * 							&&, elements in the inner sets by ||)
 * @return	set of time coordinates where an operator can be placed to start and end on a timeslice
 * 			boundary, and fulfill the required conditions given by the parameters
 */
std::set<t_coord> find_possible_operator_time_coords(t_extent operator_t_extent,
		const timeslice_sizes& timeslice_sizes_pattern, const std::set<std::set<int> >& no_boundary_times) {

	std::set<t_coord> possible_times;

	for (int operator_pos_t = 0; operator_pos_t < accumulate(timeslice_sizes_pattern.begin(), timeslice_sizes_pattern.end(), 0);
			operator_pos_t = next_timeslice_boundary(operator_pos_t, timeslice_sizes_pattern)) {

		bool non_boundaries_fulfilled = true;
		for (const int one_of_t_not_a_boundary : no_boundary_times) {
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

void try_place_operator(const timeslice_sizes& timeslice_sizes_pattern, const std::set<t_extent> t_extents,
		std::map<timeslice_sizes, std::map<t_extent, std::vector<t_coord> > >& possible_operator_placements) {

	std::set<std::set<int> > no_boundary_times;
	//TODO construct no_boundary_times

	for (const t_extent operator_t_extent : t_extents) {
		std::set<t_coord> possible_operator_t_coords = find_possible_operator_time_coords(operator_t_extent,
				timeslice_sizes_pattern, no_boundary_times);
		//TODO store results in possible_operator_pacements
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

	//TODO change program syntax, parse new parameters
	const int tsl_num = 2, min_size = 1, max_size = 3;
	const set<int> ts = { 4, 5, 6 };
	const vector<pair<int, int> > size_limits(tsl_num, { min_size, max_size + 1 });

	map<timeslice_sizes, map<t_extent, vector<t_coord> > > possible_operator_placements;
	nest_for(size_limits,
			try_place_operator, ts, possible_operator_placements);

	//TODO report results

/*
	vector<int> no_boundary_times;
	if (!handle_GNU_options(argc, argv, no_boundary_times))
		return 0;

	const vector<int> timeslice_sizes = parse_unsigned_int_list(argv[1]);
	const int t_extent = stoi(argv[2]);

	set<int> possible_ts = find_possible_operator_time_coords(t_extent, timeslice_sizes, no_boundary_times);
	cout << possible_ts.size() << " possible positions: ";
	for (const int t : possible_ts)
		cout << t << " ";
	cout << "\n";
	*/

}
