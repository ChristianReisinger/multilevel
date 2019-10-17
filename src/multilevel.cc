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

#include <fields.hh>
#include <io.hh>
#include <heatbath.hh>	//indirect dependency artifact ..
#include <smearing_techniques.hh>
#include <geometry2.hh>
#include <LinkPath.hh>
#include <io_tools.hh>
#include <helper_functions.hh>
#include <linear_algebra.hh>

#include <sublattice_algebra.hh>
#include <T_field.hh>
#include <MultilevelConfig.hh>
#include <MultilevelAnalyzer.hh>
#include <LevelDef.hh>
#include <twolink_operator_functions.hh>
#include <OperatorFactor.hh>
#include <TwolinkOperator.hh>
#include <TwolinkComputer.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

using latticetools_0719::SUN_elems;

void print_help(char* argv0) {
	std::cout << "Usage: " << argv0
			<< " [-m] [-e <ext>] [-b <beta> -s <seed> -u <level_updates> [-w]] <T> <L> <WL_Rs> <NAPEs> <level_config_num> <composition_file> <config_prefix> <config_id>\n"
					"\n"
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
					"\t\tGauge config filename without id extension (.01, .01.01, etc.)\n"
					"\n"
					"\t<config_id>\n"
					"\t\tID of the top-level config\n"
					"\n"
					"\t<outfile>\n"
					"\t\tOutput filename\n"
					"\n"
					"Options\n"
					"\n"
					"\t-m\n"
					"\t\tshow required memory and exit\n"
					"\n"
					"\t-e <ext>\n"
					"\t\tappend '.<ext>' to all output file names\n"
					"\n"
					"\t-b <beta> -s <seed> -u <level_updates> [-w]\n"
					"\t\tThese options must be used together. When used, configs are generated during the multilevel algorithm\n"
					"\t\tvia heatbath with <beta>, <seed> and <level_updates>. <level_updates> is a comma separated list of the\n"
					"\t\tnumber of updates at each level in order from highest to lowest level; updates at level 0 are applied\n"
					"\t\tonce to the initial config from file before computing observables.\n"
					"\t\tWhen using also -w, generated configs are written to file with '.multilevel' appended to filenames.\n";
}

void handle_GNU_options(int argc, char**& argv, bool& show_mem,
		bool& generate, bool& write, double& beta, int& seed, std::vector<int>& level_updates,
		std::string& extension) {
	static struct option long_opts[] = {
			{ "mem", no_argument, 0, 'm' },
			{ "beta", required_argument, 0, 'b' },
			{ "seed", required_argument, 0, 's' },
			{ "updates", required_argument, 0, 'u' },
			{ "write", no_argument, 0, 'w' },
			{ "extension", no_argument, 0, 'e' },
			{ 0, 0, 0, 0 }
	};

	int opt = -1, long_opts_i = 0;
	while ((opt = getopt_long(argc, argv, "mb:s:u:we:", long_opts, &long_opts_i)) != -1) {
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
			case 'w':
				generate = true;
				write = true;
			break;
			case 'e':
				extension = optarg;
			break;
		}
	}
	argv = argv + optind - 1;
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

std::vector<TwolinkComputer> make_twolink_computers()
{
	const std::vector<std::tuple<std::string, twolink_operator_sig (*), int> > twolink_comp_args = {
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

template<typename T>
std::vector<OperatorFactor*> parse_operator_factors(const std::vector<T>& available_operators, const std::string& factor_str) {
	std::vector<OperatorFactor*> factors;
	std::regex factor_format("(\\S+)");
	for (std::sregex_iterator factor_it(factor_str.begin(), factor_str.end(), factor_format);
			factor_it != std::sregex_iterator(); ++factor_it) {

		const std::string factor_name = factor_it->str(1);
		const auto& named_factor_it = std::find_if(available_operators.begin(), available_operators.end(), [&](const T& f) {
			return f.name() == factor_name;
		});

		if (named_factor_it != available_operators.end()) {
			OperatorFactor* f = &*named_factor_it;
			factors.push_back(f);
		} else
			throw std::invalid_argument("Twolink operator " + factor_name + " does not exist");
	}
	return factors;
}

std::vector<bool> parse_timeslice_defined(const std::string& timeslice_def_str) {
	std::vector<bool> timeslice_defined;
	std::istringstream is_def_iss(timeslice_def_str);
	char is_def_char;
	while (is_def_iss >> is_def_char) {
		if (is_def_char == '.')
			timeslice_defined.push_back(false);
		else if (is_def_char == 'x')
			timeslice_defined.push_back(true);
	}
	return timeslice_defined;
}

std::pair<bool, std::sregex_iterator> make_operators_it(const std::string& operators_str) {
	const std::string name_tsldef_regex = "(\\S+?):([.x]+):";
	const std::string descr_T_regex = "(.*):(\\d+):"; //TODO T is not needed anymore .. computed automatically in TwolinkOperator class
	const std::string factors_regex = "((?:(?: |\t)+\\S+)+)\n";
	std::regex operator_format(name_tsldef_regex + factors_regex);
	std::regex operator_format_top(name_tsldef_regex + descr_T_regex + factors_regex);

	bool top = false;
	if (std::regex_search(operators_str, operator_format_top))
		top = true;

	return {top, std::sregex_iterator(operators_str.begin(), operators_str.end(), top ? operator_format_top : operator_format)};
}

template<typename T>
void parse_operators(LevelDef& level, const std::vector<T>& available_operators, const std::string& operators_str) {

	auto istop_it = make_operators_it(operators_str);
	for (auto operator_it = istop_it.second; operator_it != std::sregex_iterator(); ++operator_it) {

		const std::string operator_name = operator_it->str(1);

		TwolinkOperator curr_op(operator_name,
				parse_timeslice_defined(operator_it->str(2)),
				parse_operator_factors(available_operators, operator_it->str(istop_it.first ? 5 : 3)));

		if (istop_it.first)
			curr_op.descr(operator_it->str(3));

		if (istop_it.first || std::count_if(level.operators().begin(), level.operators().end(), [](const TwolinkOperator& op) {
			return op.name() == operator_name;
		}) == 0)
			level.add_operator(curr_op);
		else
			throw std::runtime_error("duplicate definition of '" + operator_name + "'");
	}
}

bool verify_timeslice_sizes(const std::vector<LevelDef>& levels) {
	if (levels.empty())
		return false;
	for (auto level_it = ++levels.begin(); level_it != levels.end(); ++level_it) {

		const auto& tsl_sizes = level_it->timeslice_sizes();
		auto size_it = tsl_sizes.begin();
		for (const int prev_level_size : std::prev(level_it)->timeslice_sizes()) {
			int curr_total_size = 0;

			while (size_it != tsl_sizes.end()) {
				curr_total_size += *size_it;
				++size_it;
				if (curr_total_size == prev_level_size)
					break;
				else if (curr_total_size > prev_level_size)
					return false;
			}
			if (curr_total_size != prev_level_size)
				return false;
		}
		if (size_it != tsl_sizes.end())
			return false;
	}
	return true;
}

std::vector<LevelDef> parse_levels(const std::vector<TwolinkComputer>& twolink_computers, int T, const std::string& compstr) {
	std::vector<LevelDef> levels;

	std::regex format(""
			"(thickness \\d+(,\\d+)*\n+"
			"((?!thickness)\\S+?:[.x]+:(( |\t)+\\S+)+\n+)+)+"
			"thickness T\n+"
			"((?!thickness)\\S+?:[.x]+:.*:\\d+:(( |\t)+\\S+)+\n+)+", std::regex::nosubs);
	if (!std::regex_match(compstr, format))
		throw std::runtime_error("invalid composition format");

	std::regex level_format(""
			"thickness (\\d+(?:,\\d+)*|T)\n+"
			"((?:(?!thickness)\\S+?:[.x]+(?::.*:\\d+)?:(?:(?: |\t)+\\S+)+\n+)+)");

	const std::sregex_iterator level_begin(compstr.begin(), compstr.end(), level_format);
	if (std::distance(level_begin, std::sregex_iterator()) < 2)
		throw std::invalid_argument("less than 2 levels");

	for (auto level_it = level_begin; level_it != std::sregex_iterator(); ++level_it) {

		std::vector<int> timeslice_sizes;
		const std::string size_str = level_it->str(1);
		if (size_str == "T")
			timeslice_sizes.push_back(T);
		else
			timeslice_sizes = tools::helper::parse_unsigned_int_list(size_str.c_str());
		levels.insert(levels.begin(), LevelDef(timeslice_sizes));

		if (level_it == level_begin)
			parse_operators(levels[0], twolink_computers, level_it->str(2));
		else
			parse_operators(levels[0], levels[1].operators(), level_it->str(2));
	}
	if (!verify_timeslice_sizes(levels))
		throw std::invalid_argument("invalid timeslice sizes");

	return levels;
}

void open_outfiles(std::map<std::string, std::unique_ptr<std::ofstream> >& outfiles,
		const std::map<std::string, std::pair<std::string, std::string> >& operator_filename_lineprefix) {

	for (const auto& e : operator_filename_lineprefix) {
		const std::string& filename = e.second.first;

		if (!outfiles.count(filename)) {
			if (tools::io_tools::file_exists(filename))
				throw std::runtime_error("output file '" + filename + "' already exists");

			outfiles[filename] = tools::helper::make_unique<std::ofstream>(filename);
			if (outfiles.at(filename)->fail())
				throw std::runtime_error("could not open output file '" + filename + "'");

			*outfiles.at(filename) << std::scientific << std::setprecision(11) << std::setfill(' ');
		}
	}
}

}
}
}

int main(int argc, char** argv) {
	using namespace std;
	using de_uni_frankfurt_itp::reisinger::latticetools_0719::LinkPath;
	using de_uni_frankfurt_itp::reisinger::tools::helper::parse_unsigned_int_list;
	using de_uni_frankfurt_itp::reisinger::tools::io_tools::file_exists;
	using de_uni_frankfurt_itp::reisinger::tools::helper::make_unique;
	using namespace de_uni_frankfurt_itp::reisinger::multilevel_0819;

	auto start_time = chrono::steady_clock::now();

	if (argc < 8) {
		print_help(argv[0]);
		return 1;
	}

	const auto twolink_computers = make_twolink_computers();

//	Parameters ****************************************************************************************************************************

	int T, L, config_lv0_id, seed = 0;
	double beta = 0.0;
	bool show_mem = false, generate = false, write = false;
	set<int> WL_Rs, NAPEs;
	string outfile_extension = "";
	vector<LevelDef> levels;

	try {
		vector<int> level_updates;
		handle_GNU_options(argc, argv, show_mem, generate, write, beta, seed, level_updates, outfile_extension);
		if (generate && (beta <= 0.0 || seed <= 0))
			throw std::invalid_argument("invalid <beta> or <seed>");

		T = std::stoi(argv[1]);
		L = std::stoi(argv[2]);

		vector<int> WL_R_list = parse_unsigned_int_list(argv[3]);
		WL_Rs = set<int>(WL_R_list.begin(), WL_R_list.end());

		vector<int> NAPE_list = parse_unsigned_int_list(argv[4]);
		NAPEs = set<int>(NAPE_list.begin(), NAPE_list.end());

		vector<int> level_config_num = parse_unsigned_int_list(argv[5]);

		ifstream compositions_ifs(argv[6]);
		ostringstream compstr_oss;
		compstr_oss << compositions_ifs.rdbuf();
		levels = parse_levels(twolink_computers, T, compstr_oss.str());

		if (level_config_num.size() != levels.size()
				|| (generate && level_updates.size() != levels.size()))
			throw std::invalid_argument("invalid <level_config_num> or <level_updates>");
		for (int lv_i = 0; lv_i <= levels.size(); ++lv_i) {
			levels[lv_i].config_num(level_config_num[lv_i]);
			if (generate)
				levels[lv_i].update_num(level_updates[lv_i]);
		}

		config_lv0_id = std::stoi(argv[8]);
	} catch (std::exception& e) {
		cerr << "Error: " + e.what() + "\n";
		return 1;
	}

	if (show_mem) {
		cout << "This computation uses "
				<< memory_used(levels, T, L, WL_Rs.size()) / 1000000.0
				<< " MB memory.\n";
		return 0;
	}

//	***************************************************************************************************************************************

	MultilevelConfig multilevel_config(argv[7], config_lv0_id, T, L, levels, beta, seed, write);
	MultilevelAnalyzer multilevel(multilevel_config, WL_Rs, levels);

	map<string, map<int, T_field> > T_fields;
	try {
		T_fields = multilevel.compute_T_fields();
	} catch (runtime_error& e) {
		cout << "Error: " << e.what() << "\n";
		return 1;
	}

	double* gauge_field, *smeared_gauge_field;
	multilevel_config.get(gauge_field);
	Gauge_Field_Alloc_silent(&smeared_gauge_field, T, L);
	Gauge_Field_Copy(smeared_gauge_field, gauge_field, T, L);

//	***************************************************************************************************************************************

	map<string, map<string, vector<pair<int, complex> > > > results;

	int smeared_steps = 0;
	for (int NAPE : NAPEs) {
		for (; smeared_steps < NAPE; ++smeared_steps)
			APE_Smearing_Step(smeared_gauge_field, T, L, 0.5);

		for (const auto& name_rfields : T_fields) {
			for (const auto& r_field : name_rfields.second) {
				const int WL_R = r_field.first;
				const string& fieldname = name_rfields.first;

				complex WL_avg;
				co_eq_zero(&WL_avg);

				const auto& ts = r_field.second.defined_ts();
				for (int t : ts) {
					for (int x = 0; x < L; ++x) {
						for (int y = 0; y < L; ++y) {
							for (int z = 0; z < L; ++z) {
								for (int i = 1; i < 4; ++i) {
									LinkPath S0(smeared_gauge_field, T, L, { t + operator_T_Toffset[fieldname].second, x, y, z });
									LinkPath ST(smeared_gauge_field, T, L,
											{ t + operator_T_Toffset[fieldname].first + operator_T_Toffset[fieldname].second,
													x, y, z });
									for (int r = 0; r < WL_R; ++r) {
										S0(i, true);
										ST(i, true);
									}
									complex curr_WL;
									close_Wilson_loop(&curr_WL, r_field.second.T_at(t, x, y, z, i), S0.path, ST.path);
									co_pl_eq_co(&WL_avg, &curr_WL);
								}
							}
						}
					}
				}

				co_di_eq_re(&WL_avg, 3.0 * L * L * L);
				const string filename = operator_filename_lineprefix.at(fieldname).first;
				std::ostringstream params_oss;
				params_oss << NAPE << " " << WL_R << " " << operator_filename_lineprefix.at(fieldname).second;
				results[filename][params_oss.str()].emplace_back(ts.size(), WL_avg);
			}
		}
	}

//	***************************************************************************************************************************************

	map<string, std::unique_ptr<ofstream> > outfiles;
	try {
		open_outfiles(outfiles, operator_filename_lineprefix);
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

	Gauge_Field_Free(&smeared_gauge_field);

	cout << "\nComputation time\n"
			"\tfull program : " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start_time).count() << " s\n"
			"\tgenerating configs : " << multilevel_config.milliseconds_spent_generating() / 1000 << " s\n"
			"\tcomputing observables : " << multilevel.milliseconds_spent_computing() / 1000 << " s\n";

	return 0;
}
