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
#include <exception>
#include <regex>
#include <memory>

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
#include <twolink_operators.hh>
#include <MultilevelConfig.hh>
#include <MultilevelAnalyzer.hh>

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
					"\t\t<line_prefix> (can be empty) is written at the start of the line of the corresponding output file.\n"
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

void handle_GNU_options(int argc, char**& argv,
		bool& show_mem, bool& generate, bool& write, double& beta, int& seed, std::vector<int>& level_updates,
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

double memory_used(const std::vector<std::map<std::string, std::vector<bool> > >& level_operator_timeslice_defined,
		int T, int L, int num_R) {
	const double gauge_field_bytes = T * L * L * L * 4 * SUN_elems * sizeof(double);
	const double bytes_per_timeslice = L * L * L * 3 * SO_elems * sizeof(double);

	double timeslice_num = 0.0;
	for (const auto& level : level_operator_timeslice_defined)
		for (const auto& operator_tsl : level)
			for (bool is_defined : operator_tsl.second)
				if (is_defined)
					timeslice_num += 1.0;

	return num_R * timeslice_num * bytes_per_timeslice + 3 * gauge_field_bytes;
}

void parse_operator_factors(std::vector<std::string>& operator_def,
		const std::set<std::string>& prev_level_operator_names, const std::string& factor_str) {

	std::regex factor_format("(\\S+)");
	for (std::sregex_iterator factor_it(factor_str.begin(), factor_str.end(), factor_format);
			factor_it != std::sregex_iterator(); ++factor_it) {

		const std::string factor_name = factor_it->str(1);
		if (prev_level_operator_names.count(factor_name))
			operator_def.push_back(factor_name);
		else
			throw std::runtime_error("invalid operator name '" + factor_name + "'");
	}
}

void parse_timeslice_defined(std::vector<bool>& timeslice_defined, const std::string& timeslice_def_str) {
	std::istringstream is_def_iss(timeslice_def_str);
	char is_def_char;
	while (is_def_iss >> is_def_char) {
		if (is_def_char == '.')
			timeslice_defined.push_back(false);
		else if (is_def_char == 'x')
			timeslice_defined.push_back(true);
	}
}

void parse_operators(
		std::map<std::string, std::vector<std::string> >& operator_factors,
		std::map<std::string, std::vector<bool> >& operator_timeslice_defined,
		std::map<std::string, std::pair<int, int> >& operator_T_Toffset,
		std::map<std::string, std::pair<std::string, std::string> >& operator_filename_lineprefix,
		const std::string& outfile_extension,
		const bool top,
		const std::set<std::string>& prev_level_operator_names,
		const std::string& operators_str) {

	std::regex operator_format(std::string("(\\S+?):([.x]+):")
			+ (top ? "(.*):(\\d+):" : "")
			+ "((?:(?: |\t)+\\S+)+)\n");
	for (std::sregex_iterator operator_it(operators_str.begin(), operators_str.end(), operator_format);
			operator_it != std::sregex_iterator(); ++operator_it) {

		std::vector<std::string> operator_def;
		const std::string factors = operator_it->str(top ? 5 : 3);
		parse_operator_factors(operator_def, prev_level_operator_names, factors);

		std::vector<bool> timeslice_defined;
		const std::string timeslice_def_str = operator_it->str(2);
		parse_timeslice_defined(timeslice_defined, timeslice_def_str);

		const std::string operator_name = operator_it->str(1);
		if (operator_factors.count(operator_name))
			throw std::runtime_error("duplicate definition of '" + operator_name + "'");

		operator_factors[operator_name] = operator_def;
		operator_timeslice_defined[operator_name] = timeslice_defined;
		if (top) {
			operator_filename_lineprefix[operator_name] = {
				operator_name + (outfile_extension.empty() ? "": ".") + outfile_extension,
				operator_it->str(3)};
			operator_T_Toffset[operator_name] = {tools::io_tools::parse_int(operator_it->str(4)), 0};
		}
	}
}

bool verify_thickness(const std::vector<int>& timeslice_thickness, const std::vector<int>& prev_level_thickness, bool top) {
	for (int th : timeslice_thickness)
		if (th < 1)
			return false;

	if (!prev_level_thickness.empty()) {
		auto curr_level_it = timeslice_thickness.begin();
		auto prev_level_it = prev_level_thickness.begin();
		int prev_lev_total = 0;
		while (prev_level_it != prev_level_thickness.end()) {
			prev_lev_total += *(prev_level_it++);
			if (prev_lev_total == *curr_level_it) {
				prev_lev_total = 0;
				++curr_level_it;
			} else if (prev_lev_total > *curr_level_it)
				return false;
		}
		if (top)
			return (*curr_level_it % prev_lev_total == 0);
		else
			return (curr_level_it == timeslice_thickness.end());
	}

	return true;
}

void parse_compositions(
		std::vector<std::vector<int> >& level_thickness,
		std::vector<std::map<std::string, std::vector<std::string> > >& level_operator_factors,
		std::vector<std::map<std::string, std::vector<bool> > >& level_operator_timeslice_defined,
		std::map<std::string, std::pair<int, int> >& operator_T_Toffset,
		std::map<std::string, std::pair<std::string, std::string> >& operator_filename_lineprefix,
		std::set<std::string> prev_level_operator_names,
		const std::string& outfile_extension, int T, const std::string& compstr) {

	level_thickness.clear();
	level_operator_factors.clear();
	level_operator_timeslice_defined.clear();
	operator_T_Toffset.clear();
	operator_filename_lineprefix.clear();

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
	for (std::sregex_iterator level_it(compstr.begin(), compstr.end(), level_format);
			level_it != std::sregex_iterator(); ++level_it) {

		bool top = false;
		std::vector<int> timeslice_thickness;
		const std::string thickness_val = level_it->str(1);
		if (thickness_val == "T") {
			timeslice_thickness.push_back(T);
			top = true;
		} else
			timeslice_thickness = tools::helper::parse_unsigned_int_list(thickness_val.c_str());

		if (!verify_thickness(timeslice_thickness,
				level_thickness.empty() ? std::vector<int>() : level_thickness.at(0), top))
			throw std::runtime_error("invalid timeslice partition");

		level_thickness.insert(level_thickness.begin(), timeslice_thickness);

		std::map<std::string, std::vector<std::string> > curr_operator_factors;
		std::map<std::string, std::vector<bool> > curr_operator_timeslice_defined;
		parse_operators(curr_operator_factors, curr_operator_timeslice_defined,
				operator_T_Toffset, operator_filename_lineprefix, outfile_extension,
				top, prev_level_operator_names, level_it->str(2));
		level_operator_factors.insert(level_operator_factors.begin(), curr_operator_factors);
		level_operator_timeslice_defined.insert(level_operator_timeslice_defined.begin(), curr_operator_timeslice_defined);

		prev_level_operator_names.clear();
		for (const auto& name_factors : curr_operator_factors)
			prev_level_operator_names.insert(name_factors.first);
	}
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

	if (argc < 7) {
		print_help(argv[0]);
		return 1;
	}
	const vector<string> lowest_level_operator_names = {
			"U_x_U", "UU_x_UU",
			"Ex_x_I", "ex_x_I", "Ey_x_I", "ey_x_I", "Ez_x_I", "ez_x_I",
			"Bx_x_I", "bx_x_I", "By_x_I", "by_x_I", "Bz_x_I", "bz_x_I",
			"I_x_Ex", "I_x_ex", "I_x_Ey", "I_x_ey", "I_x_Ez", "I_x_ez",
			"I_x_Bx", "I_x_bx", "I_x_By", "I_x_by", "I_x_Bz", "I_x_bz"
	};

	int T, L, config_lv0_id, seed = 0;
	double beta = 0.0;
	bool show_mem = false, generate = false, write = false;
	set<int> WL_Rs, NAPEs;

	vector<int> level_config_num, level_updates;
	vector<vector<int> > level_thickness;
	vector<map<string, vector<string> > > level_operator_factors;
	vector<map<string, vector<bool> > > level_operator_timeslice_defined;
	map<string, pair<int, int> > operator_T_Toffset;

	map<string, pair<string, string> > operator_filename_lineprefix;
	string outfile_extension = "";

	handle_GNU_options(argc, argv, show_mem, generate, write, beta, seed, level_updates, outfile_extension);
	if (generate) {
		if (beta <= 0.0) {
			cerr << "Error: invalid <beta>\n";
			return 1;
		}
		if (seed <= 0) {
			cerr << "Error: invalid <seed>\n";
			return 1;
		}
		if (level_updates.size() == 0) {
			cerr << "Error: invalid <level_updates>\n";
			return 1;
		}
	}

	stringstream arg_ss;
	arg_ss << argv[1] << " " << argv[2] << " " << argv[8];
	arg_ss >> T >> L >> config_lv0_id;
	vector<int> WL_R_list = parse_unsigned_int_list(argv[3]);
	WL_Rs = set<int>(WL_R_list.begin(), WL_R_list.end());

	vector<int> NAPE_list = parse_unsigned_int_list(argv[4]);
	NAPEs = set<int>(NAPE_list.begin(), NAPE_list.end());

	stringstream config_id_ss;
	if (arg_ss.fail()) {
		cerr << "Error: failed to read one or more of <T> <L> <config_id>\n";
		return 1;
	}

	level_config_num = parse_unsigned_int_list(argv[5]);
	if (level_config_num.size() < 2) {
		cerr << "Error: less than 2 levels\n";
		return 1;
	}

	ifstream compositions_ifs(argv[6]);
	ostringstream compstr_oss;
	compstr_oss << compositions_ifs.rdbuf();
	try {
		parse_compositions(level_thickness, level_operator_factors, level_operator_timeslice_defined,
				operator_T_Toffset, operator_filename_lineprefix,
				set<string>(lowest_level_operator_names.begin(), lowest_level_operator_names.end()),
				outfile_extension, T, compstr_oss.str());
	} catch (const std::exception& e) {
		cerr << "Error reading composition file: " << e.what() << "\n";
		return 1;
	}

	if (level_config_num.size() != level_thickness.size()) {
		cerr << "Error: invalid <level_config_num>\n";
		return 1;
	}

	if (generate && level_updates.size() != level_thickness.size()) {
		cerr << "Error: invalid <level_updates>\n";
		return 1;
	}

	if (level_operator_factors.size() != level_thickness.size()) {
		cerr << "Error: invalid compositions\n";
		return 1;
	}

	if (show_mem) {
		cout << "This computation uses "
				<< memory_used(level_operator_timeslice_defined, T, L, WL_Rs.size()) / 1000000.0
				<< " MB memory.\n";
		return 0;
	}

//	***************************************************************************************************************************************

	map<string, void (*)(double*, const double*, int, int, int&, int, int, int, int, int)> lowest_level_functions;
	int i = 0;
	for (auto& f : {
			U_x_U, UU_x_UU,
			Ex_x_I, Exbar_x_I, Ey_x_I, Eybar_x_I, Ez_x_I, Ezbar_x_I,
			Bx_x_I, Bxbar_x_I, By_x_I, Bybar_x_I, Bz_x_I, Bzbar_x_I,
			I_x_Ex, I_x_Exbar, I_x_Ey, I_x_Eybar, I_x_Ez, I_x_Ezbar,
			I_x_Bx, I_x_Bxbar, I_x_By, I_x_Bybar, I_x_Bz, I_x_Bzbar }) {
		lowest_level_functions[lowest_level_operator_names.at(i)] = f;
		++i;
	}

	MultilevelConfig multilevel_config(argv[7], config_lv0_id, T, L, level_thickness, level_config_num, beta, seed, level_updates, write);
	MultilevelAnalyzer multilevel(multilevel_config, WL_Rs, level_operator_factors, level_operator_timeslice_defined,
			lowest_level_functions);

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

	map<string, std::unique_ptr<ofstream> > outfiles;
	try {
		open_outfiles(outfiles, operator_filename_lineprefix);
	} catch (runtime_error& e) {
		cerr << "Error: " << e.what() << "\n";
		return 1;
	}

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

				int timeslice_num = 0;
				for (int t : r_field.second.defined_ts()) {
					++timeslice_num;
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

				co_di_eq_re(&WL_avg, 3.0 * L * L * L * timeslice_num);
				const string filename = operator_filename_lineprefix.at(fieldname).first;
				const string lineprefix = operator_filename_lineprefix.at(fieldname).second;
				*outfiles.at(filename) << NAPE << " " << WL_R << " " << lineprefix << " "
						<< showpos << WL_avg.re << " " << WL_avg.im << noshowpos << "\n";
			}
		}
	}

	Gauge_Field_Free(&smeared_gauge_field);

	cout << "\nComputation time\n"
			"\tfull program : " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start_time).count() << " s\n"
			"\tgenerating configs : " << multilevel_config.milliseconds_spent_generating() / 1000 << " s\n"
			"\tcomputing observables : " << multilevel.milliseconds_spent_computing() / 1000 << " s\n";

	return 0;
}
