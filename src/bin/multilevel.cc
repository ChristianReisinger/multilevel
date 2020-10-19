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
#include <Stopwatch.hh>

#include <algorithm/FactorInterface.hh>
#include <algorithm/LevelDef.hh>
#include <algorithm/MultilevelAnalyzer.hh>
#include <algorithm/MultilevelConfig.hh>
#include <algorithm/T_field.hh>
#include <algorithm/TwolinkComputer.hh>
#include <algorithm/TwolinkOperator.hh>
#include <log/logger.hh>
#include <parameters/ConfigParameters.hh>
#include <parameters/parse_parameters.hh>
#include <physics/sublattice_algebra.hh>
#include <physics/twolink_operator_functions.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

using latticetools_0719::SUN_elems;

void print_syntax_help(const char* argv0) {
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
		"\t\t\t(sizes <sz>\n"
		"\t\t\t(<name>:<timeslice_def>: <factor> ...\n"
		"\t\t\t)...\n"
		"\t\t\t)...\n"
		"\t\t\tsizes T\n"
		"\t\t\t(<file>:<timeslice_def>:<line_prefix>:<T>: <factors> ...\n"
		"\t\t\t)...\n"
		"\t\twhere levels are defined in order from lowest to highest and the top level is indicated by <sz> = 'T'\n"
		"\t\t('sizes T').\n"
		"\t\t<sz> is a comma separated list of timeslice sizes. If the total sizes does not fill the\n"
		"\t\t\tlattice, the pattern is repeated periodically. Each level must define a partition of each timeslice\n"
		"\t\t\tdefined at the next-lowest level.\n"
		"\t\t<timeslice_def> is a string of '.' or 'x', with length equal to the number of entries in <sz>, where the\n"
		"\t\t\tn-th character 'x' / '.' indicates that the operator is / is not defined at the n-th timeslice. At top\n"
		"\t\tlevel (level 0), <timeslice_def> corresponds to entries in <sz> at level 1 instead.\n"
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
		"\t\tshow approximate required memory and exit. Omit the parameters <config_prefix>, <config_id> when using\n"
		"\t\tthis option.\n"
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
		"\t\tWhen using also -w, generated configs are written to file with '.multilevel' appended to filenames.\n"
		"\n";
}

bool handle_GNU_options(int& argc, char**& argv,
		Settings& settings, ConfigParameters& config_params, std::vector<int>& level_updates) {

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
				settings.show_mem = true;
				break;
			case 'b':
				settings.generate = true;
				config_params.beta = std::stod(optarg);
				break;
			case 's':
				settings.generate = true;
				config_params.seed = std::stoi(optarg);
				break;
			case 'u':
				settings.generate = true;
				level_updates = tools::helper::parse_unsigned_int_list(optarg);
				break;
			case 'r':
				settings.generate = true;
				config_params.overrelax_steps = std::stoi(optarg);
				break;
			case 'w':
				settings.generate = true;
				config_params.write = true;
				break;
			case 'e':
				settings.outfile_extension = optarg;
				break;
			case 'h':
				print_syntax_help(argv[0]);
				print_option_help();
				return false;
				break;
		}
	}

	if (settings.generate && (config_params.beta <= 0.0 || config_params.seed <= 1))
		throw std::invalid_argument("invalid <beta> or <seed>");

	argc -= optind;
	argv = argv + optind - 1;

	return true;
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

bool handle_required_args(int argc, const char*const* argv, const std::vector<int>& level_updates,
		Settings& settings, ConfigParameters& config_params, std::vector<LevelDef>& levels) {
	using tools::helper::parse_unsigned_int_list;

	const int REQUIRED_ARG_NUM = 8 - (settings.show_mem ? 2 : 0);
	if (argc != REQUIRED_ARG_NUM) {
		print_syntax_help(argv[0]);
		return false;
	}

	const auto twolink_computers = make_twolink_computers();

	config_params.T = std::stoi(argv[1]);
	config_params.L = std::stoi(argv[2]);

	const std::vector<int> WL_R_list = parse_unsigned_int_list(argv[3]);
	settings.WL_Rs = std::set<int>(WL_R_list.begin(), WL_R_list.end());

	std::vector<int> NAPE_list = parse_unsigned_int_list(argv[4]);
	settings.NAPEs = std::set<int>(NAPE_list.begin(), NAPE_list.end());

	const std::vector<int> level_config_num = parse_unsigned_int_list(argv[5]);
	if (level_config_num.at(0) != 1)
		throw std::invalid_argument("number of global configs != 1 not implemented");

	std::ifstream compositions_ifs(argv[6]);
	std::ostringstream compstr_oss;
	compstr_oss << compositions_ifs.rdbuf();
	levels = parse_parameters::levels(twolink_computers, compstr_oss.str(), config_params.T);

	if (level_config_num.size() != levels.size()
		|| (settings.generate && level_updates.size() != levels.size()))
		throw std::invalid_argument("invalid <level_config_num> or <level_updates>");

	for (size_t lv_i = 0; lv_i < levels.size(); ++lv_i) {
		levels[lv_i].config_num(level_config_num[lv_i]);
		if (settings.generate)
			levels[lv_i].update_num(level_updates[lv_i]);
	}

	if (!settings.show_mem) {
		config_params.filestem = argv[7];
		config_params.config_lv0_id = std::stoi(argv[8]);
	}

	return true;
}

struct WLTimesliceAverage {
	std::size_t measured_timeslice_num;
	complex value;
};
using ParameterString = std::string;
// a WL can be distributed differently over timeslices (e.g. '| U U | U |' or '| U | U U |')
// -> collect WLTimesliceAverage's for equivalent WL's and take average later
using WLResults = std::map<ParameterString, std::vector<WLTimesliceAverage> >;
using Filename = std::string;
using OutputData = std::map<Filename, WLResults>;

OutputData compute_wilson_loops(const std::vector<LevelDef>& levels, const MultilevelConfig& config,
		const ConfigParameters& config_params, const Settings& settings) {

	using namespace latticetools_0719;
	using tools::Stopwatch;

	OutputData results;

	double* smeared_gauge_field;
	Gauge_Field_Alloc(smeared_gauge_field, config_params.T, config_params.L);
	Gauge_Field_Copy(smeared_gauge_field, config.get(), config_params.T, config_params.L);

	int smeared_steps = 0;
	for (int NAPE : settings.NAPEs) {
		Stopwatch APE_watch;
		std::cout << logger::timestamp() << "(IIIa) Smearing spatial links, " << NAPE << " steps ... " << std::endl;
		for (; smeared_steps < NAPE; ++smeared_steps)
			APE_Smearing_Step(smeared_gauge_field, config_params.T, config_params.L, 0.5);
		std::cout << logger::timestamp() << "(IIIa) finished in " << APE_watch.check().count() << " ms\n";

		for (const auto& op : levels[0].operators()) {
			for (const int WL_R : settings.WL_Rs) {
				complex WL_avg;
				co_eq_zero(&WL_avg);

				const auto& ts = op.defined_ts(WL_R);
				for (int t : ts) {
					complex WLs_at_t[config_params.spatial_volume() * 3];
#pragma omp parallel for collapse(3)
					for (int x = 0; x < config_params.L; ++x) {
						for (int y = 0; y < config_params.L; ++y) {
							for (int z = 0; z < config_params.L; ++z) {
								for (int i = 1; i < 4; ++i) {
									LinkPath S0(smeared_gauge_field, config_params.T, config_params.L, { t, x, y, z });
									LinkPath ST(smeared_gauge_field, config_params.T, config_params.L, { t + op.t_extent(), x, y, z });
									for (int r = 0; r < WL_R; ++r) {
										S0(i, true);
										ST(i, true);
									}
									complex curr_WL;
									double T_op[SO_elems];
									op.at(T_op, t, x, y, z, i, WL_R);
									close_Wilson_loop(&curr_WL, T_op, S0.path, ST.path);
									WLs_at_t[i - 1 + 3 * (z + config_params.L * (y + config_params.L * x))] = curr_WL; //avoid critical region
								}
							}
						}
					}
					for (int i = 0; i < config_params.spatial_volume() * 3; ++i)
						co_pl_eq_co(&WL_avg, &WLs_at_t[i]);
				}

				co_di_eq_re(&WL_avg, 3.0 * config_params.spatial_volume());
				const std::string filename = op.name() + "." + settings.outfile_extension;
				const std::string params = std::to_string(NAPE) + " " + std::to_string(WL_R) + " " + op.descr();
				results[filename][params].push_back( { ts.size(), WL_avg });
			}
		}
	}
	Gauge_Field_Free(smeared_gauge_field);

	return results;
}

auto open_outfiles(const std::set<Filename>& filenames) {
	std::map<Filename, std::unique_ptr<std::ofstream> > outfiles;
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
	return outfiles;
}

void write_results(const OutputData& results) {
	std::set<Filename> filenames;
	for (const auto& result : results)
		filenames.insert(result.first);

	auto outfiles = open_outfiles(filenames);

	for (const auto& filename_WLresults : results) {
		const auto& filename = filename_WLresults.first;
		for (const auto& params_tslaverages : filename_WLresults.second) {
			const auto& params = params_tslaverages.first;
			complex avg { 0, 0 };
			int tsl_num = 0;
			for (const auto& tsl_avg : params_tslaverages.second) {
				tsl_num += tsl_avg.measured_timeslice_num;
				co_pl_eq_co(&avg, &tsl_avg.value);
			}
			co_di_eq_re(&avg, tsl_num);
			*outfiles.at(filename) << params << " " << std::showpos << avg.re << " " << avg.im << std::noshowpos << "\n";
		}
	}
}

}
}
}

int main(int argc, char** argv) {
	using namespace std;
	using namespace de_uni_frankfurt_itp::reisinger;
	using namespace multilevel_0819;
	using tools::io_tools::file_exists;
	using tools::helper::make_unique;
	using tools::Stopwatch;

	for (int i = 0; i < argc; ++i)
		cout << argv[i] << " ";
	cout << "\n\n";

	Stopwatch program_watch;

//	Options & arguments *******************************************************************************************************************

	ConfigParameters config_params;
	Settings settings;

	vector<LevelDef> levels;

	try {
		vector<int> level_updates;
		if (!handle_GNU_options(argc, argv, settings, config_params, level_updates)) return 0;
		if (!handle_required_args(argc, argv, level_updates, settings, config_params, levels)) return 0;
	} catch (std::exception& e) {
		cerr << "Error: " << e.what() << "\n";
		return 1;
	}

	if (settings.show_mem) {
		cout << "This computation uses "
			<< memory_used(levels, config_params.T, config_params.L, settings.WL_Rs.size()) / 1000000.0
			<< " MB memory.\n";
		return 0;
	}

//	***************************************************************************************************************************************

#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
			cout << "Using " << omp_get_num_threads() << " OMP threads\n";
	}

	Stopwatch intermediate_watch;
	cout << logger::timestamp() << "(I) Initializing multilevel algorithm ... \n";
	MultilevelConfig multilevel_config(config_params);
	cout << logger::timestamp() << "(I) finished in " << intermediate_watch.reset().count() << " ms\n";

//	***************************************************************************************************************************************

	cout << logger::timestamp() << "(II) Computing temporal transporters ... " << std::endl;
	MultilevelAnalyzer multilevel(levels, multilevel_config, settings.WL_Rs);
	try {
		multilevel.compute_T_fields();
	} catch (runtime_error& e) {
		cerr << "Error: " << e.what() << "\n";
		return 1;
	}
	cout << logger::timestamp() << "(II) finished in " << intermediate_watch.reset().count() << " ms\n";

	cout << logger::timestamp() << "(III) Computing Wilson loops ... " << std::endl;
	auto results = compute_wilson_loops(levels, multilevel_config, config_params, settings);
	cout << logger::timestamp() << "(III) finished in " << intermediate_watch.reset().count() << " ms\n";

//	***************************************************************************************************************************************

	cout << logger::timestamp() << "(IV) Writing results to file ... " << std::endl;
	try {
		write_results(results);
	} catch (runtime_error& e) {
		cerr << "Error: " << e.what() << "\n";
		return 1;
	}
	cout << logger::timestamp() << "(IV) finished in " << intermediate_watch.reset().count() << " ms\n";

	cout << logger::timestamp() << "\nComputation time\n"
		"\tfull program : " << program_watch.check().count() << " ms\n"
		"\tupdating configs : " << multilevel_config.milliseconds_spent_updating() << " ms\n"
		"\tcomputing observables : " << multilevel.milliseconds_spent_computing() << " ms\n";

	return 0;
}
