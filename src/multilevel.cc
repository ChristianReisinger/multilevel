#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <getopt.h>
#include <fstream>

#include <fields.hh>
#include <io.hh>
#include <heatbath.hh>	//indirect dependency artifact ..

#include <geometry2.hh>
#include <LinkPath.hh>

#include <helper_functions.hh>

#include <linear_algebra.hh>
#include <sublattice_algebra.hh>
#include <sublattice_fields.hh>
#include <twolink_operators.hh>
#include <MultilevelAnalyzer.hh>

bool show_mem = false;

void handle_GNU_options(int argc, char**& argv, bool& generate, double& beta, int& seed) {
	static struct option long_opts[] = {
			{ "mem", no_argument, 0, 'm' },
			{ "beta", required_argument, 0, 'b' },
			{ "seed", required_argument, 0, 's' },
			{ 0, 0, 0, 0 }
	};

	int opt = -1, long_opts_i = 0;
	while ((opt = getopt_long(argc, argv, "mb:s:", long_opts, &long_opts_i)) != -1) {
		switch (opt) {
			case 'm':
				show_mem = true;
			break;
			case 'b':
				generate = true;
				std::istringstream beta_iss(optarg);
				beta_iss >> beta;
			break;
			case 's':
				generate = true;
				std::istringstream seed_iss(optarg);
				seed_iss >> seed;
			break;
		}
	}
	argv = argv + optind - 1;
}

double memory_bytes_used(std::vector<std::vector<std::vector<int> > > field_compositions, std::vector<int> level_thickness, int T,
		int L) {
	int levels = level_thickness.size();
	if (levels < 2 || field_compositions.size() != levels)
		return -1;
	double bytes_per_field = T * L * L * L * 3 * SUN_N * SUN_N * SUN_N * SUN_N * 2 * 8;
	double operators_per_site = 0.0;
	for (int level = 0; level < levels; ++level)
		operators_per_site += (double) field_compositions[level].size() / (double) level_thickness[level == 0 ? 1 : level];
	return operators_per_site * bytes_per_field;
}

int main(int argc, char **argv) {
	using namespace std;

	if (argc < 7) {
		cerr << "Usage: " << argv[0]
				<< " [-m] [-b <beta> -s <seed>] <T> <L> <level_config_num> <level_thickness> <config_prefix> <config_id> <outfile>\n";
		return 0;
	}

	bool generate = false;
	int T, L, config_lv0_id, seed = 0;
	double beta = 0.0;
	vector<int> level_thickness;
	level_thickness.push_back(T);

	handle_GNU_options(argc, argv, generate, beta, seed);
	if (generate) {
		if (beta <= 0.0) {
			cerr << "Error: invalid beta\n";
			return 0;
		}
		if (seed <= 0) {
			cerr << "Error: invalid seed\n";
			return 0;
		}
	}

	stringstream arg_ss;
	arg_ss << argv[1] << " " << argv[2] << " " << argv[6];
	arg_ss >> T >> L >> config_lv0_id;

	stringstream config_id_ss;
	if (arg_ss.fail()) {
		cerr << "Error: failed to read one or more of <T> <L> <config_id>\n";
		return 0;
	}

	vector<int> level_config_num = parse_unsigned_int_list(argv[3]);
	if (level_config_num.size() < 2) {
		cerr << "Error: less than 2 levels\n";
		return 0;
	}

	vector<int> ds = parse_unsigned_int_list(argv[4]);
	if (ds.empty()) {
		cerr << "Error: no timeslice thickness specified\n";
		return 0;
	}
	for (int d : ds)
		level_thickness.push_back(d);

//	********************************** Params **********************************

	int WL_R = 10, WL_T = 8;
	int T_offset = 0; // T_field computed with multilevel at site n corresponds to a temporal Wilson line direct product at site n + T_offset * unit_vec_t

	vector<vector<vector<int> > > field_compositions = {
			{ { 0, 0 } },
			{ { 0, 0 } }
	};

	MultilevelAnalyzer multilevel(T, L, WL_R, level_thickness, argv[5], level_config_num,
			{ UU_x_UU },
			field_compositions, generate, beta, seed, { 10, 10 });

//	****************************************************************************

	if (show_mem) {
		cout << "This computation uses " << memory_bytes_used(field_compositions, level_thickness, T, L) / 1000000.0 << " MB memory.\n";
		return 0;
	}

	ofstream out_ofs(argv[7]);
	if (out_ofs.fail()) {
		cerr << "Error: could not open output file '" << argv[7] << "'";
		return 0;
	}

//	***************************************************************************************************************************************

	const int timeslice_num = level_thickness[0] / level_thickness[1];
	double* T_field[1];
	T_field_alloc_zero(*T_field, 3, timeslice_num, L);

	multilevel.compute_sublattice_fields( { config_lv0_id }, 0, T_field);

	double* gauge_field;
	Gauge_Field_Alloc_silent(&gauge_field, T, L);
	read_gauge_field(gauge_field, multilevel.config_filename( { config_lv0_id }).c_str(), T, L);

	complex WL_avg;
	co_eq_zero(&WL_avg);
	for (int t = 0; t < T; t += level_thickness[1]) {
		for (int x = 0; x < L; ++x) {
			for (int y = 0; y < L; ++y) {
				for (int z = 0; z < L; ++z) {
					for (int i = 1; i < 4; ++i) {
						LinkPath S0(gauge_field, T, L, { t + T_offset, x, y, z }), ST(gauge_field, T, L, { t + T_offset + WL_T, x, y, z });
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
