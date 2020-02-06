#include <iostream>
#include <getopt.h>
#include <vector>
#include <set>
#include <numeric>

#include <helper_functions.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {
namespace try_operator_on_timeslices {

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

int next_timeslice_boundary(int t, const std::vector<int>& timeslice_sizes_pattern) {
	const int tsl_num = timeslice_sizes_pattern.size();
	int curr_t = 0;

	for (int tsl_i = 0; curr_t <= t; ++tsl_i)
		curr_t += timeslice_sizes_pattern[tsl_i % tsl_num];
	return curr_t;
}

bool is_timeslice_boundary(const int t, const std::vector<int>& timeslice_sizes_pattern) {
	int curr_t = 0;
	while (curr_t <= t) {
		if (t == curr_t)
			return true;

		curr_t = next_timeslice_boundary(curr_t, timeslice_sizes_pattern);
	}
	return false;
}

std::set<int> possible_operator_pos_times(int operator_t_extent,
		const std::vector<int>& timeslice_sizes, const std::vector<int>& no_boundary_times) {

	std::set<int> possible_times;

	for (int operator_pos_t = 0; operator_pos_t < accumulate(timeslice_sizes.begin(), timeslice_sizes.end(), 0);
			operator_pos_t = next_timeslice_boundary(operator_pos_t, timeslice_sizes)) {

		bool t_is_possible = true;
		for (const int t_no_boundary : no_boundary_times)
			if (is_timeslice_boundary(operator_pos_t + t_no_boundary, timeslice_sizes)) {
				t_is_possible = false;
				break;
			}

		if (t_is_possible && is_timeslice_boundary(operator_pos_t + operator_t_extent, timeslice_sizes))
			possible_times.insert(operator_pos_t);
	}

	return possible_times;
}

}
}
}
}

int main(int argc, char **argv) {
	using namespace std;
	using namespace de_uni_frankfurt_itp::reisinger;
	using namespace multilevel_0819::try_operator_on_timeslices;
	using tools::helper::parse_unsigned_int_list;

	vector<int> no_boundary_times;
	if (!handle_GNU_options(argc, argv, no_boundary_times))
		return 0;

	const vector<int> timeslice_sizes = parse_unsigned_int_list(argv[1]);
	const int t_extent = stoi(argv[2]);

	set<int> possible_ts = possible_operator_pos_times(t_extent, timeslice_sizes, no_boundary_times);
	cout << possible_ts.size() << " possible positions";
	if (!possible_ts.empty()) {
		cout << " - times ";
		for (const int t : possible_ts)
			cout << t << " ";
	}
	cout << "\n";

}
