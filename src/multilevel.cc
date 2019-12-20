#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <getopt.h>
#include <chrono>
#include <stdexcept>
#include <regex>
#include <memory>
#include <numeric>
#include <algorithm>
#include <omp.h>

#include <smearing_techniques.hh>
#include <linear_algebra.hh>
#include <geometry2.hh>
#include <LinkPath.hh>
#include <io_tools.hh>
#include <helper_functions.hh>

#include <sublattice_algebra.hh>
#include <T_field.hh>
#include <MultilevelConfig.hh>
#include <MultilevelAnalyzer.hh>
#include <LevelDef.hh>
#include <twolink_operator_functions.hh>
#include <FactorInterface.hh>
#include <TwolinkOperator.hh>
#include <TwolinkComputer.hh>
#include <parse_parameters.hh>
#include <logger.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

using latticetools_0719::SUN_elems;

void print_syntax_help(char* argv0) {
	std::cout << "\nUsage: " << argv0 << ""
			"\t[--help] [--mem] [--extension <ext>]\n"
			"\t[--beta <beta> --seed <seed> --updates <level_updates> [--overrelax <overrelax_steps>] [--write]]\n"
			"\t<T> <L> <WL_Rs> <NAPEs> <level_config_num> <composition_file> <config_prefix> <config_id>\n\n";
}

void print_option_help() {
	std::cout << "\n"
			"Parameters\n"
			"\n"
			"\t<T> <L>\n"
			"\t\tlattice dimensions\n"
			"\n"
			"\t<WL_Rs>\n"
			"\t\tspatial Wilson loop sizes\n"
			"\n"
			"\t<NAPEs>\n"
			"\t\tlist of number of APE smearing steps. Spatial Wilson lines are smeared after computing temporal lines on\n"
			"\t\tunsmeared gauge field with multilevel\n"
			"\n"
			"\t<level_config_num>\n"
			"\t\tComma separated list of the number of configs at each level\n"
			"\n"
			"\t<composition_file>\n"
			"\t\tFile containing operator composition definitions at each level, format:\n"
			"\t\t\t(thickness <th>\n"
			"\t\t\t(<name>:<timeslice_def>: <factor> ...\n"
			"\t\t\t)...\n"
			"\t\t\t)...\n"
			"\t\t\tthickness T\n"
			"\t\t\t(<file>:<timeslice_def>:<line_prefix>:<T>: <factors> ...\n"
			"\t\t\t)...\n"
			"\t\twhere levels are defined in order from lowest to highest and the top level is indicated by <th> = 'T'\n"
			"\t\t('thickness T').\n"
			"\t\t<th> is a comma separated list of timeslice thicknesses. If the total thickness does not fill the\n"
			"\t\t\tlattice, the pattern is repeated periodically. Each level must define a partition of each timeslice\n"
			"\t\t\tdefined at the next-lowest level.\n"
			"\t\t<timeslice_def> is a string of '.' or 'x', with length equal to the number of entries in <th>, where the\n"
			"\t\t\tn-th character 'x' / '.' indicates that the operator is / is not defined at the n-th timeslice. At top\n"
			"\t\tlevel (level 0), <timeslice_def> corresponds to entries in <th> at level 1 instead.\n"
			"\t\tOperators named <name> are defined as product of <factor>'s on the same line, where <factor> is one of\n"
			"\t\t\tthe <name>'s defined at the next-lowest level.\n"
			"\t\tAt top level, '<file>.<ext>' (or '<file>' if <ext> is empty) is the output file into which final results\n"
			"\t\t\tare written.\n"
			"\t\t<T> is the temporal size of the Wilson loop for each operator.\n"
			"\t\t<line_prefix> (can be empty): output files have columns '<NAPE> <WL_R> <line_prefix> <re> <im>'.\n"
			"\t\tResults for entries with equal <file> and <line_prefix> are averaged over.\n"
			"\t\tPossible <name>'s at the lowest level are\n"
			"\t\t\tU_x_U, UU_x_UU, Ex_x_I, ex_x_I, Ey_x_I, ..., Bx_x_I, ..., I_x_Ex, ... .\n"
			"\n"
			"\t<config_prefix>\n"
			"\t\tGauge config filename without id extension (.1, .1.1, etc.)\n"
			"\n"
			"\t<config_id>\n"
			"\t\tID of the top-level config\n"
			"\n"
			"\t<outfile>\n"
			"\t\tOutput filename\n"
			"\n"
			"Options\n"
			"\n"
			"\t--memory | -m\n"
			"\t\tshow required memory and exit\n"
			"\n"
			"\t--extension -e <ext>\n"
			"\t\tappend '.<ext>' to all output file names\n"
			"\n"
			"\t--beta <beta>, -b <beta>\n"
			"\t--seed <seed>, -s <seed>\n"
			"\t--updates <level_updates>, -u <level_updates>\n"
			"\t[--write, -w]\n"
			"\t[--overrelax, -r <overrelax_steps>]\n"
			"\t\tThese options must be used together. When used, configs are generated during the multilevel algorithm\n"
			"\t\tvia heatbath with <beta>, <seed> and <level_updates>. <level_updates> is a comma separated list of the\n"
			"\t\tnumber of updates at each level in order from highest to lowest level; updates at level 0 are applied\n"
			"\t\tonce to the initial config from file before computing observables.\n"
			"\t\tSet number of overrelaxation steps with -r (default is 0), has no effect for SU(2) build.\n"
			"\t\tWhen using also -w, generated configs are written to file with '.multilevel' appended to filenames.\n";
}

bool handle_GNU_options(int argc, char**& argv, bool& show_mem,
		bool& generate, bool& write, double& beta, int& seed, std::vector<int>& level_updates, int& overrelax_steps,
		std::string& extension) {

	bool no_help_required = true;

	static struct option long_opts[] = {
			{ "mem", no_argument, 0, 'm' },
			{ "beta", required_argument, 0, 'b' },
			{ "seed", required_argument, 0, 's' },
			{ "updates", required_argument, 0, 'u' },
			{ "overrelax", required_argument, 0, 'r' },
			{ "write", no_argument, 0, 'w' },
			{ "extension", no_argument, 0, 'e' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
	};

	int opt = -1, long_opts_i = 0;
	while ((opt = getopt_long(argc, argv, "mb:s:u:r:we:h", long_opts, &long_opts_i)) != -1) {
		switch (opt) {
			case 'm':
				show_mem = true;
			break;
			case 'b':
				generate = true;
				beta = std::stod(optarg);
			break;
			case 's':
				generate = true;
				seed = std::stoi(optarg);
			break;
			case 'u':
				generate = true;
				level_updates = tools::helper::parse_unsigned_int_list(optarg);
			break;
			case 'r':
				generate = true;
				overrelax_steps = std::stoi(optarg);
			break;
			case 'w':
				generate = true;
				write = true;
			break;
			case 'e':
				extension = optarg;
			break;
			case 'h':
				print_syntax_help(argv[0]);
				print_option_help();
				no_help_required = false;
			break;
		}
	}
	argv = argv + optind - 1;

	if (argc - optind + 1 != 8) {
		print_syntax_help(argv[0]);
		no_help_required = false;
	}

	return no_help_required;
}

double memory_used(const std::vector<LevelDef>& levels, int T, int L, int num_R) {
	const double gauge_field_bytes = T * L * L * L * 4 * SUN_elems * sizeof(double);
	const double bytes_per_timeslice = L * L * L * 3 * SO_elems * sizeof(double);

	double timeslice_num = 0.0;
	for (const auto& level : levels) {

		int curr_level_tsl_num = 0;
		for (const auto& op : level.operators())
			curr_level_tsl_num += op.timeslice_num_per_cycle();

		curr_level_tsl_num *= ((double) T) / std::accumulate(level.timeslice_sizes().begin(), level.timeslice_sizes().end(), 0);
		timeslice_num += curr_level_tsl_num;
	}
	return num_R * timeslice_num * bytes_per_timeslice + 3 * gauge_field_bytes;
}

std::vector<TwolinkComputer> make_twolink_computers() {
	static const std::vector<std::tuple<std::string, twolink_operator_sig (*), int> > twolink_comp_args = {
			{ "U_x_U", U_x_U, 1 }, { "UU_x_UU", UU_x_UU, 2 },
			{ "Ex_x_I", Ex_x_I, 0 }, { "ex_x_I", Exbar_x_I, 0 },
			{ "Ey_x_I", Ey_x_I, 0 }, { "ey_x_I", Eybar_x_I, 0 },
			{ "Ez_x_I", Ez_x_I, 0 }, { "ez_x_I", Ezbar_x_I, 0 },
			{ "Bx_x_I", Bx_x_I, 0 }, { "bx_x_I", Bxbar_x_I, 0 },
			{ "By_x_I", By_x_I, 0 }, { "by_x_I", Bybar_x_I, 0 },
			{ "Bz_x_I", Bz_x_I, 0 }, { "bz_x_I", Bzbar_x_I, 0 },
			{ "I_x_Ex", I_x_Ex, 0 }, { "I_x_ex", I_x_Exbar, 0 },
			{ "I_x_Ey", I_x_Ey, 0 }, { "I_x_ey", I_x_Eybar, 0 },
			{ "I_x_Ez", I_x_Ez, 0 }, { "I_x_ez", I_x_Ezbar, 0 },
			{ "I_x_Bx", I_x_Bx, 0 }, { "I_x_bx", I_x_Bxbar, 0 },
			{ "I_x_By", I_x_By, 0 }, { "I_x_by", I_x_Bybar, 0 },
			{ "I_x_Bz", I_x_Bz, 0 }, { "I_x_bz", I_x_Bzbar, 0 }
	};
	std::vector<TwolinkComputer> twolink_comps;
	for (const auto& args : twolink_comp_args)
		twolink_comps.emplace_back(std::get<0>(args), std::get<1>(args), std::get<2>(args));
	return twolink_comps;
}

void open_outfiles(std::map<std::string, std::unique_ptr<std::ofstream> >& outfiles, const std::set<std::string>& filenames) {
	for (const auto& filename : filenames) {

		if (!outfiles.count(filename)) {
			if (tools::io_tools::file_exists(filename))
				throw std::runtime_error("output file '" + filename + "' already exists");

			outfiles[filename] = tools::helper::make_unique<std::ofstream>(filename);
			if (outfiles.at(filename)->fail())
				throw std::runtime_error("could not open output file '" + filename + "'");

			*outfiles.at(filename) << std::scientific << std::setprecision(11);
		}
	}
}

}
}
}

int main(int argc, char** argv) {
	using namespace std;
	using namespace de_uni_frankfurt_itp::reisinger;
	using namespace latticetools_0719;
	using namespace multilevel_0819;
	using tools::helper::parse_unsigned_int_list;
	using tools::io_tools::file_exists;
	using tools::helper::make_unique;

	auto start_time = chrono::steady_clock::now();
	
	const auto twolink_computers = make_twolink_computers();

//	Parameters ****************************************************************************************************************************

	int T, L, config_lv0_id, seed = 1, overrelax_steps = 0;
	double beta = 0.0;
	bool show_mem = false, generate = false, write = false;
	set<int> WL_Rs, NAPEs;
	string outfile_extension = "";
	vector<LevelDef> levels;

	try {
		vector<int> level_updates;
		if (!handle_GNU_options(argc, argv, show_mem, generate, write, beta, seed, level_updates, overrelax_steps, outfile_extension))
			return 0;

		if (generate && (beta <= 0.0 || seed <= 1))
			throw invalid_argument("invalid <beta> or <seed>");

		T = stoi(argv[1]);
		L = stoi(argv[2]);

		const vector<int> WL_R_list = parse_unsigned_int_list(argv[3]);
		WL_Rs = std::set<int>(WL_R_list.begin(), WL_R_list.end());

		vector<int> NAPE_list = parse_unsigned_int_list(argv[4]);
		NAPEs = std::set<int>(NAPE_list.begin(), NAPE_list.end());

		const vector<int> level_config_num = parse_unsigned_int_list(argv[5]);

		ifstream compositions_ifs(argv[6]);
		ostringstream compstr_oss;
		compstr_oss << compositions_ifs.rdbuf();
		levels = parse_parameters::levels(twolink_computers, compstr_oss.str(), T);

		if (level_config_num.size() != levels.size()
				|| (generate && level_updates.size() != levels.size()))
			throw invalid_argument("invalid <level_config_num> or <level_updates>");
		if (level_config_num.at(0) != 1)
			throw invalid_argument("number of configs at level 0 must be 1");
		for (size_t lv_i = 0; lv_i < levels.size(); ++lv_i) {
			levels[lv_i].config_num(level_config_num[lv_i]);
			if (generate)
				levels[lv_i].update_num(level_updates[lv_i]);
		}

		config_lv0_id = std::stoi(argv[8]);
	} catch (std::exception& e) {
		cerr << "Error: " << e.what() << "\n";
		return 1;
	}

	if (show_mem) {
		cout << "This computation uses "
				<< memory_used(levels, T, L, WL_Rs.size()) / 1000000.0
				<< " MB memory.\n";
		return 0;
	}

//	***************************************************************************************************************************************

	auto intermediate_time = chrono::steady_clock::now();
	logger::print_timestamp();
	cout << "(I) Initializing multilevel algorithm ... \n";

	MultilevelConfig multilevel_config(argv[7], config_lv0_id, T, L, beta, seed, overrelax_steps, write);
	MultilevelAnalyzer multilevel(levels, multilevel_config, WL_Rs);

	logger::print_timestamp();
	cout << "(I) finished in " << logger::get_ms_since_and_reset(intermediate_time) << " ms\n";

	logger::print_timestamp();
	cout << "(II) Computing temporal transporters ... " << std::endl;
	try {
		multilevel.compute_T_fields();
	} catch (runtime_error& e) {
		cerr << "Error: " << e.what() << "\n";
		return 1;
	}
	logger::print_timestamp();
	cout << "(II) finished in " << logger::get_ms_since_and_reset(intermediate_time) << " ms\n";

	logger::print_timestamp();
	cout << "(III) Computing Wilson loops ... " << std::endl;
	double* smeared_gauge_field;
	Gauge_Field_Alloc(smeared_gauge_field, T, L);
	Gauge_Field_Copy(smeared_gauge_field, multilevel_config.get(), T, L);

//	***************************************************************************************************************************************

	map<string, map<string, vector<pair<int, complex> > > > results;

	int smeared_steps = 0;
	for (int NAPE : NAPEs) {
		for (; smeared_steps < NAPE; ++smeared_steps)
			APE_Smearing_Step(smeared_gauge_field, T, L, 0.5);

		for (const auto& op : levels[0].operators()) {
			for (const int WL_R : WL_Rs) {
				complex WL_avg;
				co_eq_zero(&WL_avg);

				const auto& ts = op.defined_ts(WL_R);
				for (int t : ts) {
#pragma omp parallel for collapse(3)
					for (int x = 0; x < L; ++x) {
						for (int y = 0; y < L; ++y) {
							for (int z = 0; z < L; ++z) {
								for (int i = 1; i < 4; ++i) {
									LinkPath S0(smeared_gauge_field, T, L, { t, x, y, z });
									LinkPath ST(smeared_gauge_field, T, L, { t + op.t_extent(), x, y, z });
									for (int r = 0; r < WL_R; ++r) {
										S0(i, true);
										ST(i, true);
									}
									complex curr_WL;
									double T_op[SO_elems];
									op.at(T_op, t, x, y, z, i, WL_R);
									close_Wilson_loop(&curr_WL, T_op, S0.path, ST.path);
#pragma omp critical
									co_pl_eq_co(&WL_avg, &curr_WL);
								}
							}
						}
					}
				}

				co_di_eq_re(&WL_avg, 3.0 * L * L * L);
				const string filename = op.name() + "." + outfile_extension;
				const string params = to_string(NAPE) + " " + to_string(WL_R) + " " + op.descr();
				results[filename][params].emplace_back(ts.size(), WL_avg);
			}
		}
	}
	logger::print_timestamp();
	cout << "(III) finished in " << logger::get_ms_since_and_reset(intermediate_time) << " ms\n";

//	***************************************************************************************************************************************

	logger::print_timestamp();
	cout << "(IV) Writing results to file ... " << std::endl;
	map<string, unique_ptr<ofstream> > outfiles;
	try {
		set<string> filenames;
		for (const auto& result : results)
			filenames.insert(result.first);
		open_outfiles(outfiles, filenames);
	} catch (runtime_error& e) {
		cerr << "Error: " << e.what() << "\n";
		return 1;
	}

	for (const auto& filename_params_data : results) {
		const auto& filename = filename_params_data.first;
		for (const auto& params_data : filename_params_data.second) {
			const auto& params = params_data.first;
			complex avg { 0, 0 };
			int tsl_num = 0;
			for (const auto& tnum_data : params_data.second) {
				tsl_num += tnum_data.first;
				co_pl_eq_co(&avg, &tnum_data.second);
			}
			co_di_eq_re(&avg, tsl_num);
			*outfiles.at(filename) << params << " " << showpos << avg.re << " " << avg.im << noshowpos << "\n";
		}
	}

	logger::print_timestamp();
	cout << "(IV) finished in " << logger::get_ms_since_and_reset(intermediate_time) << " ms\n";

	Gauge_Field_Free(smeared_gauge_field);

	logger::print_timestamp();
	cout << "\nComputation time\n"
			"\tfull program : " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_time).count() << " ms\n"
			"\tupdating configs : " << multilevel_config.milliseconds_spent_updating() << " ms\n"
			"\tcomputing observables : " << multilevel.milliseconds_spent_computing() << " ms\n";

	return 0;
}
