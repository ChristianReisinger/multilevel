#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <getopt.h>
#include <fstream>

#include <linear_algebra.hh>
#include <geometry.hh>
#include <fields.hh>
#include <io.hh>

#include <geometry2.hh>
#include <global_defs.hh>
#include <LinkPath.hh>

#include <helper_functions.hh>

#include <sublattice_algebra.hh>
#include <sublattice_fields.hh>
#include <twolink_operators.hh>
#include <MultilevelAnalyzer.hh>

bool insert_E = false;

void handle_GNU_options(int argc, char**& argv) {
	static struct option long_opts[] = {
			{ "Ez", no_argument, 0, 'E' },
			{ 0, 0, 0, 0 }
	};

	int opt = -1, long_opts_i = 0;
	while ((opt = getopt_long(argc, argv, "E", long_opts, &long_opts_i)) != -1) {
		switch (opt) {
			case 'E':
				insert_E = true;
			break;
		}
	}
	argv = argv + optind - 1;
}

int main(int argc, char **argv) {
	using namespace std;
	int WL_R = 5, WL_T = 8;

	if (argc != 7) {
		cerr << "Usage: " << argv[0] << " <T> <L> <level_config_num> <config_prefix> <config_id> <outfile>\n";
		return 0;
	}

	int T, L, config_lv0_id;

	stringstream arg_ss;
	arg_ss << argv[1] << " " << argv[2] << " " << argv[5];
	arg_ss >> T >> L >> config_lv0_id;

	std::stringstream config_id_ss;
	if (arg_ss.fail()) {
		cerr << "Error: failed to read one or more of <T> <L> <config_id>\n";
		return 0;
	}

	vector<int> level_config_num = parse_unsigned_int_list(argv[3]);
	if (level_config_num.size() != 3) {
		cerr << "Error: invalid number of levels, 3 required\n";
		return 0;
	}

	ofstream out_ofs(argv[6]);
	if (out_ofs.fail()) {
		cerr << "Error: could not open output file '" << argv[6] << "'";
		return 0;
	}

	const vector<int> level_thickness { T, 4, 2 };
	const int timeslice_num = level_thickness[0] / level_thickness[1];

	double* T_field[1];
	T_field_alloc_zero(*T_field, 3, timeslice_num, L);

	MultilevelAnalyzer multilevel(T, L, WL_R, level_thickness, argv[4], level_config_num,
			{ compute_T_ti_T, compute_T_ti_Tclov_lower_half, compute_Tclov_upper_half_ti_T },
			{
					{ { 0, 1 } },
					{ { 0, 1 }, { 2, 0 } },
					{ { 0 }, { 1 }, { 2 } }
			});
	multilevel.compute_sublattice_fields( { config_lv0_id }, 0, T_field);

	double* gauge_field;
	Gauge_Field_Alloc(&gauge_field, T, L);
	read_gauge_field(gauge_field, multilevel.config_filename( { config_lv0_id }).c_str(), T, L);

	complex WL_avg;
	co_eq_zero(&WL_avg);
	for (int t = 0; t < T; t += level_thickness[1]) {
		for (int x = 0; x < L; ++x) {
			for (int y = 0; y < L; ++y) {
				for (int z = 0; z < L; ++z) {
					for (int i = 1; i < 4; ++i) {
						LinkPath S0(gauge_field, T, L, { t, x, y, z }), ST(gauge_field, T, L, { t + WL_T, x, y, z });
						for (int r = 0; r < WL_R; ++r) {
							S0(i, true);
							ST(i, true);
						}
						complex curr_WL;
						close_Wilson_loop(&curr_WL, *T_field + T_field_index(t, x, y, z, 3, i - 1, T, L, level_thickness[1]),
								S0.path, ST.path);
						co_pl_eq_co(&WL_avg, &curr_WL);
					}
				}
			}
		}
	}
	delete[] *T_field;
	co_di_eq_re(&WL_avg, 3.0 * L * L * L * timeslice_num);

	out_ofs << scientific << setprecision(11) << showpos << WL_avg.re << " " << WL_avg.im << "\n";
}
