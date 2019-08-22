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

void handle_GNU_options(int argc, char**& argv, bool& generate, double& beta, int& seed, std::vector<int>& level_updates) {
	static struct option long_opts[] = {
			{ "mem", no_argument, 0, 'm' },
			{ "beta", required_argument, 0, 'b' },
			{ "seed", required_argument, 0, 's' },
			{ "updates", required_argument, 0, 'u' },
			{ 0, 0, 0, 0 }
	};

	int opt = -1, long_opts_i = 0;
	while ((opt = getopt_long(argc, argv, "mb:s:u:", long_opts, &long_opts_i)) != -1) {
		std::stringstream optarg_iss;
		switch (opt) {
			case 'm':
				show_mem = true;
			break;
			case 'b':
				generate = true;
				optarg_iss << optarg;
				optarg_iss >> beta;
			break;
			case 's':
				generate = true;
				optarg_iss << optarg;
				optarg_iss >> seed;
			break;
			case 'u':
				generate = true;
				level_updates = parse_unsigned_int_list(optarg);
				level_updates.insert(level_updates.begin(), 20);
			break;
		}
	}
	argv = argv + optind - 1;
}

double memory_used(const std::vector<std::vector<std::vector<int> > >& field_compositions, const std::vector<int>& level_thickness,
		int T, int L) {
	int levels = level_thickness.size();
	double volume = T * L * L * L;
	double bytes_per_field = volume * 3 * SO_elems * sizeof(double);
	double operators_per_site = 0.0;
	for (int level = 0; level < levels; ++level)
		operators_per_site += (double) field_compositions[level].size() / (double) level_thickness[level == 0 ? 1 : level];
	double gauge_field_bytes = volume * 4 * SUN_elems * sizeof(double);
	return operators_per_site * bytes_per_field + gauge_field_bytes;
}

int main(int argc, char** argv) {
	using namespace std;

	if (argc < 7) {
		cerr << "Usage: " << argv[0]
				<< " [-m] [-b <beta> -s <seed> -u <level_updates>] <T> <L> <level_config_num> <level_thickness> <config_prefix> <config_id> <outfile>\n"
						"\n"
						"Parameters\n"
						"\n"
						"\t<T> <L>\n"
						"\t\tlattice dimensions\n"
						"\n"
						"\t<level_config_num>\n"
						"\t\tcomma separated list of the number of configs at each level\n"
						"\n"
						"\t<level_thickness>\n"
						"\t\tcomma separated list of the timeslice thickness at each level except the top one in order from highest to\n"
						"\t\tlowest level\n"
						"\n"
						"\t<config_prefix>\n"
						"\t\tgauge config filename without id extension (.01, .01.01, etc.)\n"
						"\n"
						"\t<config_id>\n"
						"\t\tid of the top-level config\n"
						"\n"
						"\t<outfile>\n"
						"\t\toutput filename\n"
						"\n"
						"Options\n"
						"\n"
						"\t-b <beta> -s <seed> -u <level_updates>\n"
						"\t\tThese options must be used together. When used, configs are generated during the multilevel algorithm\n"
						"\t\tvia heatbath with <beta>, <seed> and <level_updates>. <level_updates> is a comma separated list of the\n"
						"\t\tnumber of updates at each level except the top one in order from highest to lowest level.\n";
		return 0;
	}

	bool generate = false;
	vector<int> level_updates;
	int T, L, config_lv0_id, seed = 0;
	double beta = 0.0;

	handle_GNU_options(argc, argv, generate, beta, seed, level_updates);
	if (generate) {
		if (beta <= 0.0) {
			cerr << "Error: invalid <beta>\n";
			return 0;
		}
		if (seed <= 0) {
			cerr << "Error: invalid <seed>\n";
			return 0;
		}
		if (level_updates.size() == 0) {
			cerr << "Error: invalid <level_updates>\n";
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

	vector<int> level_thickness = parse_unsigned_int_list(argv[4]);
	level_thickness.insert(level_thickness.begin(), T);
	if (level_thickness.size() != level_config_num.size()) {
		cerr << "Error: invalid <level_thickness>\n";
		return 0;
	}

	if (generate && level_updates.size() != level_config_num.size()) {
		cerr << "Error: invalid <level_updates>\n";
		return 0;
	}

//	********************************** Params **********************************

	int WL_R = 10;

	vector<vector<vector<int> > > field_compositions = {
			/********** levels { 2 } **********/
			{ { 0, 2, 3 }, { 0, 2, 1 }, { 1, 2, 1 }, { 1, 2, 1, 3 }, { 0, 1, 2, 1, 3 }, { 0, 1, 2, 1, 2 }, { 1, 1, 2, 1, 1 } },
			/********** levels { 4, 2 } **********/
//			{ { 0, 1 }, { 0, 2 }, { 3, 2 }, { 3, 4 }, { 5, 6, 1 }, { 5, 6, 2 }, { 7, 6, 2 } },
//			{ { 0, 2 }, { 3 }, { 1 }, { 1, 2 }, { 1, 3 }, { 0, 1 }, { 2, 1 }, { 1, 1 } },
			/********** levels { 6, 2 } **********/
//			{ { 0 }, { 1 }, { 2 }, { 2, 3 }, { 4, 5 }, { 4, 6 }, { 7, 6 } },
//			{ { 0, 2, 3 }, { 0, 2, 1 }, { 1, 2, 1 }, { 3 }, { 0, 1, 2 }, { 1, 3 }, { 1, 1 }, { 1, 1, 2 } },
			/********** levels { X, 2 } **********/
			{ { 0 }, { 1 }, { 2 }, { 3 } }
	/********** *************** **********/
	/********** levels { 6, 3 } **********/
//			{ { 0 }, { 1 }, { 2, 3 }, { 2, 4 }, { 5, 4 }, { 5, 6 }, { 7, 8 } },
//			{ { 0, 1 }, { 0, 2 }, { 3, 4 }, { 1 }, { 2 }, { 5, 4 }, { 5 }, { 5, 0 }, { 5, 1 } },
//			{ { 3, 2 }, { 3 }, { 1 }, { 0, 3 }, { 2, 3 }, { 1, 3 } }
			};

	//T_Toffset[i] defines T & Toffset for field_composition[0][i]
	//T_field at site {t, x, y, z} computed with multilevel corresponds to an operator in temporal direction defined at site {t + Toffset, x, y, z}
	vector<pair<int, int> > T_Toffset = {
			/********** levels { X, 2 } **********/
			{ 4, 1 }, { 5, 1 }, { 6, 0 }, { 7, 0 }, { 8, 1 }, { 9, 1 }, { 10, 0 }
	/********** levels { 6, 3 } **********/
//			{ 4, 0 }, { 5, 0 }, { 6, 1 }, { 7, 1 }, { 8, 0 }, { 9, 0 }, { 10, 0 }
			};

	MultilevelAnalyzer multilevel(T, L, WL_R, level_thickness, argv[5], level_config_num,
			{ IU_x_IU, UU_x_UU, UU_x_UCU, U_x_U },
			field_compositions, generate, beta, seed, level_updates);

//	****************************************************************************

	if (field_compositions.size() != level_config_num.size()) {
		cerr << "Error: invalid compositions\n";
		return 0;
	}

	if (show_mem) {
		cout << "This computation uses " << memory_used(field_compositions, level_thickness, T, L) / 1000000.0 << " MB memory.\n";
		return 0;
	}

	ofstream out_ofs(argv[7]);
	if (out_ofs.fail()) {
		cerr << "Error: could not open output file '" << argv[7] << "'";
		return 0;
	}
	out_ofs << scientific << setprecision(11) << setfill(' ');

//	***************************************************************************************************************************************

	const int timeslice_num = level_thickness[0] / level_thickness[1];
	const int top_level_field_num = field_compositions[0].size();
	double* T_field[top_level_field_num];
	for (int i = 0; i < top_level_field_num; ++i)
		T_field_alloc_zero(T_field[i], 3, timeslice_num, L);

	multilevel.compute_sublattice_fields( { config_lv0_id }, 0, T_field);

	double* gauge_field;
	Gauge_Field_Alloc_silent(&gauge_field, T, L);
	read_gauge_field(gauge_field, multilevel.config_filename( { config_lv0_id }).c_str(), T, L);

	for (int T_field_i = 0; T_field_i < top_level_field_num; ++T_field_i) {
		complex WL_avg;
		co_eq_zero(&WL_avg);
		for (int t = 0; t < T; t += level_thickness[1]) {
			for (int x = 0; x < L; ++x) {
				for (int y = 0; y < L; ++y) {
					for (int z = 0; z < L; ++z) {
						for (int i = 1; i < 4; ++i) {
							LinkPath S0(gauge_field, T, L, { t + T_Toffset[T_field_i].second, x, y, z });
							LinkPath ST(gauge_field, T, L, { t + T_Toffset[T_field_i].first + T_Toffset[T_field_i].second, x, y, z });
							for (int r = 0; r < WL_R; ++r) {
								S0(i, true);
								ST(i, true);
							}
							complex curr_WL;
							close_Wilson_loop(&curr_WL, T_field[T_field_i] + T_field_index(t, x, y, z, 3, i - 1, T, L, level_thickness[1]),
									S0.path, ST.path);
							co_pl_eq_co(&WL_avg, &curr_WL);
						}
					}
				}
			}
		}
		co_di_eq_re(&WL_avg, 3.0 * L * L * L * timeslice_num);
		out_ofs << setw(2) << T_Toffset[T_field_i].first << " " << showpos << WL_avg.re << " " << WL_avg.im << noshowpos << "\n";
	}

	for (int i = 0; i < top_level_field_num; ++i)
		delete[] T_field[i];
}
