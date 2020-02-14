#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <getopt.h>

#include <helper_functions.hh>

#include <find_sizes/main.hh>

namespace de_uni_frankfurt_itp::reisinger::multilevel_0819::find_sizes {

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

}

int main(int argc, char** argv) {
	using namespace std;
	using namespace de_uni_frankfurt_itp::reisinger;
	using namespace multilevel_0819::find_sizes;
	using tools::helper::parse_unsigned_int_list;

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

	auto timeslice_setups = find_sizes(size_limits, total_pattern_sizes, t_extents, center);
	filter(timeslice_setups, required_extents);
	print(timeslice_setups);

}
