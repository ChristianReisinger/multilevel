#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <getopt.h>

#include <linear_algebra.hh>
#include <geometry.hh>
#include <fields.hh>
#include <io.hh>

#include <geometry2.hh>
#include <global_defs.hh>
#include <LinkPath.hh>

#include <helper_functions.hh>

#include <sublattice_algebra.hh>
#include <twolink_operators.hh>

//compute [ [T(0)T(1)] [T(2)T'(3)] ] [ [T''(4)T(5)] [T(6)T(7)] ]
//with T' T'' such that E_z is inserted at t=4 as a clover

int T, L;
std::string config_prefix;
std::vector<int> level_config_num;
bool insert_E = false;

std::string config_filename(const std::vector<int>& tag) {
	std::ostringstream filename_oss;
	filename_oss << config_prefix;
	for (int i = 0; i < tag.size(); ++i)
		filename_oss << "." << std::setfill('0') << std::setw(log10(level_config_num.at(i)) + 1) << tag.at(i);
	return filename_oss.str();
}

unsigned long int T_field_index(int t, int x, int y, int z, int n, int i, int T, int L, int timeslice_thickness) {
	return (n * get_index(t / timeslice_thickness, x, y, z, T / timeslice_thickness, L) + i) * SO_elems;
}

void T_field_alloc_zero(double*& T_field, int n, int timeslice_num, int L) {
	T_field = new double[n * SO_elems * L * L * L * timeslice_num]();
}

void T_field_di_eq_re(double* T_field, double re, int n, int T, int L, int timeslice_thickness) {
	for (int i = 0; i < n * SO_elems * L * L * L * T / timeslice_thickness; ++i)
		T_field[i] /= re;
}

void obtain_sublattice_gauge_field(double*& sub_gauge_field, const std::vector<int>& tag, int T, int L) {
	Gauge_Field_Alloc(&sub_gauge_field, T, L);
	read_gauge_field(sub_gauge_field, config_filename(tag).c_str(), T, L);
}

void multilevel(const std::vector<int>& conf_tag, const std::vector<int>& level_thickness, int level, int WL_R, double** T_fields) {
	if (level == 0) {
		double* level1_T_fields[2];
		for (int i = 0; i < 2; ++i)
			T_field_alloc_zero(level1_T_fields[i], 3, T / level_thickness[1], L);
		multilevel(conf_tag, level_thickness, level + 1, WL_R, level1_T_fields);

		for (int t = 0; t < level_thickness[0]; t += level_thickness[1]) {
			for (int x = 0; x < L; ++x) {
				for (int y = 0; y < L; ++y) {
					for (int z = 0; z < L; ++z) {
						for (int i = 0; i < 3; ++i) {
							int sub_index_first = T_field_index(t, x, y, z, 3, i, T, L, level_thickness[1]);
							int sub_index_second = T_field_index(t + level_thickness[level + 1], x, y, z, 3, i, T, L, level_thickness[1]);

							so_eq_so_ti_so(*T_fields + sub_index_first,
									level1_T_fields[0] + sub_index_first,
									level1_T_fields[1] + sub_index_second);
						}
					}
				}
			}
		}
		for (int i = 0; i < 2; ++i)
			delete[] level1_T_fields[i];

	} else if (level == 1) {

		for (int conf = 1; conf <= level_config_num[level]; ++conf) {
			std::vector<int> curr_tag(conf_tag);
			curr_tag.push_back(conf);

			double* level2_T_fields[3];
			for (int i = 0; i < 3; ++i)
				T_field_alloc_zero(level2_T_fields[i], 3, T / level_thickness[2], L);
			multilevel(curr_tag, level_thickness, level + 1, WL_R, level2_T_fields);

			for (int t = 0; t < level_thickness[0]; t += level_thickness[level]) {
				for (int x = 0; x < L; ++x) {
					for (int y = 0; y < L; ++y) {
						for (int z = 0; z < L; ++z) {
							for (int i = 0; i < 3; ++i) {
								int curr_index = T_field_index(t, x, y, z, 3, i, T, L, level_thickness[level]);
								int sub_index_first = T_field_index(t, x, y, z, 3, i, T, L, level_thickness[level + 1]);
								int sub_index_second = T_field_index(t + level_thickness[level + 1], x, y, z, 3, i, T, L,
										level_thickness[level + 1]);
								double U[2][SO_elems];

								so_eq_so_ti_so(U[0], level2_T_fields[0] + sub_index_first, level2_T_fields[1] + sub_index_second);
								so_eq_so_ti_so(U[1], level2_T_fields[2] + sub_index_first, level2_T_fields[0] + sub_index_second);

								for (int j = 0; j < 2; ++j)
									so_pl_eq_so(T_fields[j] + curr_index, U[j]);
							}
						}
					}
				}
			}
			for (int i = 0; i < 3; ++i)
				delete[] level2_T_fields[i];
		}
		for (int i = 0; i < 2; ++i)
			T_field_di_eq_re(T_fields[i], level_config_num[level], 3, T, L, level_thickness[level]);

	} else if (level == 2) {
		for (int conf = 1; conf <= level_config_num[level]; ++conf) {
			std::vector<int> curr_tag(conf_tag);
			curr_tag.push_back(conf);

			double* sub_gauge_field;
			obtain_sublattice_gauge_field(sub_gauge_field, curr_tag, level_thickness[0], L);

			for (int t = 0; t < T; t += level_thickness[level]) {
				for (int x = 0; x < L; ++x) {
					for (int y = 0; y < L; ++y) {
						for (int z = 0; z < L; ++z) {
							for (int i = 1; i < 4; ++i) {
								int curr_index = T_field_index(t, x, y, z, 3, i - 1, T, L, level_thickness[level]);
								double U[3][SO_elems];
								compute_T_ti_T(U[0], sub_gauge_field, t, x, y, z, i, WL_R);
								compute_T_ti_Tclov_lower_half(U[1], sub_gauge_field, t, x, y, z, i, WL_R);
								compute_Tclov_upper_half_ti_T(U[2], sub_gauge_field, t, x, y, z, i, WL_R);

								for (int j = 0; j < 3; ++j)
									so_pl_eq_so(T_fields[j] + curr_index, U[j]);
							}
						}
					}
				}
			}
			Gauge_Field_Free(&sub_gauge_field);
		}
		for (int i = 0; i < 3; ++i)
			T_field_di_eq_re(T_fields[i], level_config_num[level], 3, T, L, level_thickness[level]);
	}

}

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

	if (argc != 6) {
		cerr << "Usage: " << argv[0] << " [-E | --Ez] <T> <L> <level_config_num> <config_prefix> <config_id>\n";
		return 0;
	}

	int config_lv0_id;

	stringstream arg_ss;
	arg_ss << argv[1] << " " << argv[2] << " " << argv[5];
	arg_ss >> T >> L >> config_lv0_id;

	std::stringstream config_id_ss;
	if (arg_ss.fail()) {
		cerr << "Error: failed to read one or more of <T> <L> <config_id>\n";
		return 0;
	}

	level_config_num = parse_unsigned_int_list(argv[3]);
	if (level_config_num.size() != 3) {
		cerr << "Error: invalid number of levels, 3 required\n";
		return 0;
	}

	config_prefix = argv[4];

	double* gauge_field;
	Gauge_Field_Alloc(&gauge_field, T, L);
	read_gauge_field(gauge_field, config_filename( { config_lv0_id }).c_str(), T, L);

	const std::vector<int> level_thickness { T, 4, 2 };
	const int timeslice_num = level_thickness[0] / level_thickness[1];

	double* T_field[1];
	T_field_alloc_zero(*T_field, 3, timeslice_num, L);
	multilevel( { config_lv0_id }, level_thickness, 0, WL_R, T_field);

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

	cout << scientific << setprecision(11) << showpos << WL_avg.re << " " << WL_avg.im << "\n";
}
