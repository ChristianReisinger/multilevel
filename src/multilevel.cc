#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <getopt.h>
#include <fstream>
#include <chrono>
#include <exception>
#include <regex>

#include <fields.hh>
#include <io.hh>
#include <heatbath.hh>	//indirect dependency artifact ..

#include <geometry2.hh>
#include <LinkPath.hh>

#include <io_tools.hh>
#include <helper_functions.hh>

#include <linear_algebra.hh>
#include <sublattice_algebra.hh>
#include <sublattice_fields.hh>
#include <twolink_operators.hh>
#include <MultilevelConfig.hh>
#include <MultilevelAnalyzer.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

using latticetools_0719::SUN_elems;

bool show_mem = false;

void print_help(char* argv0) {
	std::cout << "Usage: " << argv0
			<< " [-m] [-b <beta> -s <seed> -u <level_updates> [-w]] <T> <L> <WL_R> <level_config_num> <composition_file> <config_prefix> <config_id>\n"
					"\n"
					"Parameters\n"
					"\n"
					"\t<T> <L>\n"
					"\t\tlattice dimensions\n"
					"\n"
					"\t<WL_R>\n"
					"\t\tspatial Wilson loop size\n"
					"\n"
					"\t<level_config_num>\n"
					"\t\tComma separated list of the number of configs at each level\n"
					"\n"
					"\t<composition_file>\n"
					"\t\tFile containing operator composition definitions at each level, format:\n"
					"\t\t\t(thickness <th>\n"
					"\t\t\t(<label> <factor> ...\n"
					"\t\t\t)...\n"
					"\t\t\t)...\n"
					"\t\t\tthickness T\n"
					"\t\t\t(<filename>:<line_prefix>:<T>: <factors> ...\n"
					"\t\t\t)...\n"
					"\t\twhere levels with thickness <th> are defined in order from lowest to highest and the top level is\n"
					"\t\tindicated by 'thickness T'. Operators named <label> are defined as product of <factor>'s on the same\n"
					"\t\tline, where <factor> is on of the <label>'s from the next-lowest level. At the top level, <filename>\n"
					"\t\tis the output file into which the final results are written, <T> is the temporal size of the Wilson\n"
					"\t\tloop for each operator and <line_prefix> (can be empty) is written in the output at the start of the\n"
					"\t\tline for the corresponding result. Possible <label>'s at the lowest level are\n"
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
					"\t-b <beta> -s <seed> -u <level_updates> [-w]\n"
					"\t\tThese options must be used together. When used, configs are generated during the multilevel algorithm\n"
					"\t\tvia heatbath with <beta>, <seed> and <level_updates>. <level_updates> is a comma separated list of the\n"
					"\t\tnumber of updates at each level in order from highest to lowest level; updates at level 0 are applied\n"
					"\t\tonce to the initial config from file before computing observables.\n"
					"\t\tWhen using also -w, generated configs are written to file with '.multilevel' appended to filenames.\n";
}

void handle_GNU_options(int argc, char**& argv, bool& generate, bool& write, double& beta, int& seed, std::vector<int>& level_updates) {
	static struct option long_opts[] = {
			{ "mem", no_argument, 0, 'm' },
			{ "beta", required_argument, 0, 'b' },
			{ "seed", required_argument, 0, 's' },
			{ "updates", required_argument, 0, 'u' },
			{ "write", no_argument, 0, 'w' },
			{ 0, 0, 0, 0 }
	};

	int opt = -1, long_opts_i = 0;
	while ((opt = getopt_long(argc, argv, "mb:s:u:w", long_opts, &long_opts_i)) != -1) {
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
		}
	}
	argv = argv + optind - 1;
}

double memory_used(
		const std::vector<std::vector<std::vector<int> > >& field_compositions,
		const std::vector<int>& level_thickness,
		int T, int L) {
	const int levels = level_thickness.size();
	const double volume = T * L * L * L;
	const double bytes_per_field = volume * 3 * SO_elems * sizeof(double);
	double operators_per_site = 0.0;
	for (int level = 0; level < levels; ++level)
		operators_per_site += (double) field_compositions[level].size() / (double) level_thickness[level == 0 ? 1 : level];
	const double gauge_field_bytes = volume * 4 * SUN_elems * sizeof(double);
	return operators_per_site * bytes_per_field + 2 * gauge_field_bytes;
}

void parse_operator_factors(std::vector<int>& operator_def,
		const std::map<std::string, int>& operator_labels,
		const std::string& factor_str) {

	std::regex factor_format("(\\S+)");
	for (std::sregex_iterator factor_it(factor_str.begin(), factor_str.end(), factor_format);
			factor_it != std::sregex_iterator(); ++factor_it) {

		const std::string factor_label = factor_it->str(1);
		if (operator_labels.count(factor_label))
			operator_def.push_back(operator_labels.at(factor_label));
		else {
			std::ostringstream err_iss;
			err_iss << "invalid operator name '" << factor_label << "'";
			throw std::runtime_error(err_iss.str());
		}
	}
}

void parse_operators(
		std::vector<std::vector<int> >& operators,
		std::map<std::string, int>& labels,
		std::vector<std::pair<std::string, std::string> >& operator_filename_lineprefix,
		std::vector<std::pair<int, int> >& T_Toffset,
		const std::map<std::string, int>& prev_level_labels,
		const bool top,
		const std::string& operators_str) {

	std::ostringstream operator_format_oss;
	operator_format_oss << "(\\S+?)"
			<< (top ? ":(.*):(\\d+):" : "")
			<< "((?:(?: |\t)+\\S+)+)\n";
	std::regex operator_format(operator_format_oss.str());
	for (std::sregex_iterator operator_it(operators_str.begin(), operators_str.end(), operator_format);
			operator_it != std::sregex_iterator(); ++operator_it) {

		std::vector<int> operator_def;
		const std::string factors = operator_it->str(top ? 4 : 2);
		parse_operator_factors(operator_def, prev_level_labels, factors);

		const std::string operator_name = operator_it->str(1);
		const int operator_index = operators.size();
		labels[operator_name] = operator_index;
		operators.push_back(operator_def);
		if (top) {
			operator_filename_lineprefix.push_back( { operator_name, operator_it->str(2) });
			T_Toffset.push_back( { tools::io_tools::parse_int(operator_it->str(3)), 0 });
		}
	}
}

void parse_compositions(
		std::vector<int>& level_thickness,
		std::vector<std::vector<std::vector<int> > >& field_compositions,
		std::vector<std::pair<int, int> >& T_Toffset,
		std::vector<std::pair<std::string, std::string> >& operator_filename_lineprefix,
		const std::string& compstr, int T) {
	std::map<std::string, int> lowest_level_labels;
	int label_i = 0;
	for (std::string label : { "U_x_U", "UU_x_UU",
			"Ex_x_I", "ex_x_I", "Ey_x_I", "ey_x_I", "Ez_x_I", "ez_x_I",
			"Bx_x_I", "bx_x_I", "By_x_I", "by_x_I", "Bz_x_I", "bz_x_I",
			"I_x_Ex", "I_x_ex", "I_x_Ey", "I_x_ey", "I_x_Ez", "I_x_ez",
			"I_x_Bx", "I_x_bx", "I_x_By", "I_x_by", "I_x_Bz", "I_x_bz" }) {
		lowest_level_labels[label] = label_i;
		++label_i;
	}

	std::regex format(""
			"(thickness \\d+\n+"
			"((?!thickness)\\S+(( |\t)+\\S+)+\n+)+)+"
			"thickness T\n+"
			"((?!thickness)\\S+?:.*:\\d+:(( |\t)+\\S+)+\n+)+", std::regex::nosubs);
	if (!std::regex_match(compstr, format))
		throw std::runtime_error("invalid composition format");

	std::regex level_format(""
			"thickness (\\d+|T)\n+"
			"((?:(?!thickness)\\S+?(?::.*:\\d+:)?(?:(?: |\t)+\\S+)+\n+)+)");
	bool lowest = true;
	std::map<std::string, int> prev_level_labels;
	for (std::sregex_iterator level_it(compstr.begin(), compstr.end(), level_format);
			level_it != std::sregex_iterator(); ++level_it) {

		const std::map<std::string, int>& operator_labels = (lowest ? lowest_level_labels : prev_level_labels);

		bool top = false;
		int thickness = -1;
		const std::string thickness_val = level_it->str(1);
		if (thickness_val == "T") {
			thickness = T;
			top = true;
		} else
			thickness = tools::io_tools::parse_int(thickness_val);
		if (thickness < 1 || (!lowest &&
				(thickness <= level_thickness.at(0) || thickness % level_thickness.at(0) != 0)
				))
			throw std::runtime_error("invalid level thickness");
		level_thickness.insert(level_thickness.begin(), thickness);

		std::vector<std::vector<int> > curr_level_operators;
		std::map<std::string, int> curr_level_labels;

		parse_operators(curr_level_operators, curr_level_labels, operator_filename_lineprefix, T_Toffset,
				operator_labels, top, level_it->str(2));

		prev_level_labels = curr_level_labels;
		field_compositions.insert(field_compositions.begin(), curr_level_operators);

		lowest = false;
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
		return 0;
	}

	int T, L, WL_R, config_lv0_id, seed = 0;
	double beta = 0.0;
	bool generate = false, write = false;

	vector<int> level_config_num, level_updates, level_thickness;
	vector<vector<vector<int> > > field_compositions;
	vector<pair<int, int> > T_Toffset;
	vector<pair<string, string> > operator_filename_lineprefix;

	handle_GNU_options(argc, argv, generate, write, beta, seed, level_updates);
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
	arg_ss << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[7];
	arg_ss >> T >> L >> WL_R >> config_lv0_id;

	stringstream config_id_ss;
	if (arg_ss.fail()) {
		cerr << "Error: failed to read one or more of <T> <L> <config_id>\n";
		return 0;
	}

	level_config_num = parse_unsigned_int_list(argv[4]);
	if (level_config_num.size() < 2) {
		cerr << "Error: less than 2 levels\n";
		return 0;
	}

	ifstream compositions_ifs(argv[5]);
	ostringstream compstr_oss;
	compstr_oss << compositions_ifs.rdbuf();
	try {
		parse_compositions(level_thickness, field_compositions, T_Toffset, operator_filename_lineprefix, compstr_oss.str(), T);
	} catch (const std::exception& e) {
		cerr << "Error when reading composition file: '" << e.what() << "'\n";
		return 0;
	}

	if (level_config_num.size() != level_thickness.size()) {
		cerr << "Error: invalid <level_config_num>\n";
		return 0;
	}

	if (generate && level_updates.size() != level_thickness.size()) {
		cerr << "Error: invalid <level_updates>\n";
		return 0;
	}

	if (field_compositions.size() != level_thickness.size()) {
		cerr << "Error: invalid compositions\n";
		return 0;
	}

	if (show_mem) {
		cout << "This computation uses " << memory_used(field_compositions, level_thickness, T, L) / 1000000.0
				<< " MB memory.\n";
		return 0;
	}

	map<string, std::unique_ptr<ofstream> > outfiles;
	for (const auto& e : operator_filename_lineprefix) {
		const string filename = e.first;
		if (!outfiles.count(filename)) {
			if (file_exists(filename)) {
				cerr << "Error: output file '" << filename << "' already exists\n";
				return 0;
			}
			outfiles[filename] = make_unique<ofstream>(filename);
			if (outfiles.at(filename)->fail()) {
				cerr << "Error: could not open output file '" << filename << "'";
				return 0;
			}
			*outfiles.at(filename) << scientific << setprecision(11) << setfill(' ');
		}
	}

//	***************************************************************************************************************************************
	const vector<void (*)(double*, const double*, int, int, int&, int, int, int, int, int)> lowest_level_functions = {
			U_x_U, UU_x_UU,
			Ex_x_I, Exbar_x_I, Ey_x_I, Eybar_x_I, Ez_x_I, Ezbar_x_I,
			Bx_x_I, Bxbar_x_I, By_x_I, Bybar_x_I, Bz_x_I, Bzbar_x_I,
			I_x_Ex, I_x_Exbar, I_x_Ey, I_x_Eybar, I_x_Ez, I_x_Ezbar,
			I_x_Bx, I_x_Bxbar, I_x_By, I_x_Bybar, I_x_Bz, I_x_Bzbar
	};

	MultilevelConfig multilevel_config(argv[6], config_lv0_id, T, L, level_thickness, level_config_num, beta, seed, level_updates, write);
	MultilevelAnalyzer multilevel(multilevel_config, WL_R, field_compositions, lowest_level_functions);

	const int timeslice_num = level_thickness[0] / level_thickness[1];
	const int top_level_field_num = field_compositions[0].size();
	double* T_fields[top_level_field_num];
	for (int i = 0; i < top_level_field_num; ++i)
		T_field_alloc_zero(T_fields[i], 3, timeslice_num, L);

	multilevel.compute_sublattice_fields(T_fields);
	multilevel_config.update(0);
	double* gauge_field;
	multilevel_config.get(gauge_field);

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
							close_Wilson_loop(&curr_WL, T_fields[T_field_i] + T_field_index(t, x, y, z, 3, i - 1, T, L, level_thickness[1]),
									S0.path, ST.path);
							co_pl_eq_co(&WL_avg, &curr_WL);
						}
					}
				}
			}
		}
		co_di_eq_re(&WL_avg, 3.0 * L * L * L * timeslice_num);
		const string filename = operator_filename_lineprefix.at(T_field_i).first;
		const string lineprefix = operator_filename_lineprefix.at(T_field_i).second;
		*outfiles.at(filename) << lineprefix << " " << showpos << WL_avg.re << " " << WL_avg.im << noshowpos << "\n";
	}

	for (int i = 0; i < top_level_field_num; ++i)
		delete[] T_fields[i];

	cout << "\nComputation time\n"
			"\tfull program : " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start_time).count() << " s\n"
			"\tgenerating configs : " << multilevel_config.milliseconds_spent_generating() / 1000 << " s\n"
			"\tcomputing observables : " << multilevel.milliseconds_spent_computing() / 1000 << " s\n";
}
