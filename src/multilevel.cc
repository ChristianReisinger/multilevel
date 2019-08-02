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

//compute [ [T(0)T(1)] [T(2)T'(3)] ] [ [T''(4)T(5)] [T(6)T(7)] ]
//with T' T'' such that E_z is inserted at t=4 as a clover

int T, L;
std::string config_prefix;
std::vector<int> level_config_num { 1 };
bool insert_E = false;

std::string config_filename(const std::vector<int>& tag) {
	std::ostringstream filename_oss;
	filename_oss << config_prefix;
	for (int i : tag)
		filename_oss << "." << i;
	return filename_oss.str();
}

/**
 * n : starting lattice site of the full operator
 * t_sub : initial time of the timeslice on which to compute the sublattice operator
 */
void sublattice_operator(const std::vector<int>& conf_tag, int t_sub, const std::array<int, 4>& n, int rsep, int dir, double* ret) {

	double* sub_gauge_field;
	Gauge_Field_Alloc(&sub_gauge_field, T, L);
	read_gauge_field(sub_gauge_field, config_filename(conf_tag).c_str(), T, L);

	double T0[SUN_elems], TR[SUN_elems];
	LinkPath p(sub_gauge_field, T, L, n);
	cm_eq_cm(T0, p.move(0, t_sub)(0, true)(0, true).path);
	cm_eq_cm(TR, p.reset(n).move(0, t_sub).move(dir, rsep)(0, true)(0, true).path);

	if (insert_E) {
		double clov[SUN_elems], clov_dag[SUN_elems], U[SUN_elems];
		if (t_sub == 2) {
			cm_eq_cm(clov, p.reset(n).move(0, t_sub + 2).move(dir, rsep)(dir, true)(0, false)(dir, false)(0, true).path);
			cm_pl_eq_cm(clov, p.reset(n).move(0, t_sub + 2).move(dir, rsep)(0, false)(dir, false)(0, true)(dir, true).path);

			cm_eq_cm_dag(clov_dag, clov);
			cm_ti_eq_re(clov_dag, -1);
			cm_pl_eq_cm(clov, clov_dag);
			cm_ti_eq_re(clov, 0.5);

			cm_eq_cm_ti_cm(U, TR, clov);
			cm_eq_cm(TR, U);
		} else if (t_sub == 4) {
			cm_eq_cm(clov, p.reset(n).move(0, t_sub).move(dir, rsep)(0, true)(dir, true)(0, false)(dir, false).path);
			cm_pl_eq_cm(clov, p.reset(n).move(0, t_sub).move(dir, rsep)(dir, false)(0, true)(dir, true)(0, false).path);

			cm_eq_cm_dag(clov_dag, clov);
			cm_ti_eq_re(clov_dag, -1);
			cm_pl_eq_cm(clov, clov_dag);
			cm_ti_eq_re(clov, 0.5);

			cm_eq_cm_ti_cm(U, clov, TR);
			cm_eq_cm(TR, U);
		}
	}

	so_eq_cm_x_cm(ret, T0, TR);
}
void multilevel(const std::vector<int>& conf_tag, int WL_T, int level,
		int t_offset, const std::array<int, 4>& n, int dir, int rsep,
		double* const ret) {

	const std::vector<int> level_thickness { WL_T, 4, 2 };

	so_eq_zero(ret);
	for (int conf = 1; conf <= level_config_num[level]; ++conf) {
		std::vector<int> curr_tag(conf_tag);
		if (level != 0)
			curr_tag.push_back(conf);

		double curr_so[SO_elems];

		if (level < 2) {
			so_eq_id(curr_so);
			for (int t_sub = 0; t_sub < level_thickness[level]; t_sub += level_thickness[level + 1]) {
				double subavg[SO_elems], U[SO_elems];
				multilevel(curr_tag, WL_T, level + 1, t_offset + t_sub, n, dir, rsep, subavg);
				so_eq_so_ti_so(U, curr_so, subavg);
				so_eq_so(curr_so, U);
			}
		} else {
			sublattice_operator(curr_tag, t_offset, n, rsep, dir, curr_so);
		}

		so_pl_eq_so(ret, curr_so);
	}
	so_ti_eq_re(ret, 1.0 / level_config_num[level]);
}

/**
 * compute the Wilson loop from the sublattice operator SO and spatial Wilson lines
 * S0 at t = 0, and ST at t = T = temporal Wilson loop size
 */
complex complete_Wilson_loop(const double* SO, const double* S0, const double* ST) {
	double S0_ti_T[SUN_elems];
	cm_eq_zero(S0_ti_T);

	for (int a = 0; a < SUN_N; ++a) {
		for (int b = 0; b < SUN_N; ++b) {
			for (int i = 0; i < SUN_N; ++i) {
				for (int j = 0; j < SUN_N; ++j) {
					const double* S0_elem_ptr = S0 + cm_superindex(i, j);
					complex S0_elem { *S0_elem_ptr, *(S0_elem_ptr) };
					const double* T_elem_ptr = SO + so_superindex(i, a, j, b);
					complex T_elem { *T_elem_ptr, *(T_elem_ptr) };

					double* S0_ti_T_elem_ptr = S0_ti_T + cm_superindex(a, b);
					complex S0_ti_T_elem;
					co_eq_co_ti_co(&S0_ti_T_elem, &S0_elem, &T_elem);
					*S0_ti_T_elem_ptr += S0_ti_T_elem.re;
					*(S0_ti_T_elem_ptr + 1) += S0_ti_T_elem.im;
				}
			}
		}
	}
	double S0_ti_T_ti_ST[SUN_elems];
	cm_eq_cm_ti_cm(S0_ti_T_ti_ST, S0_ti_T, ST);

	complex WL;
	co_eq_tr_cm(&WL, S0_ti_T_ti_ST);
	co_di_eq_re(&WL, SUN_N);
	return WL;
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
	int WL_r = 5, WL_T = 8;

	if (argc != 6) {
		cerr << "Usage: " << argv[0] << " <T> <L> <level_config_num> <config_prefix> <config_id>\n";
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

	vector<int> lv1_lv2_config_num = parse_unsigned_int_list(argv[3]);
	if (lv1_lv2_config_num.size() != 2) {
		cerr << "Error: invalid number of levels, 2 required\n";
		return 0;
	}
	level_config_num.push_back(lv1_lv2_config_num[0]);
	level_config_num.push_back(lv1_lv2_config_num[1]);

	config_prefix = argv[4];

	double* gauge_field;
	Gauge_Field_Alloc(&gauge_field, T, L);
	read_gauge_field(gauge_field, config_filename( { config_lv0_id }).c_str(), T, L);

	complex WL_avg;
	co_eq_zero(&WL_avg);

	for (int t = 0; t < T; t += 2) {
		for (int x = 0; x < L; ++x) {
			for (int y = 0; y < L; ++y) {
				for (int z = 0; z < L; ++z) {
					std::array<int, 4> n = { t, x, y, z };
					for (int i = 1; i < 4; ++i) {

						double SO[SO_elems];
						multilevel( { config_lv0_id }, WL_T, 0, 0, n, i, WL_r, SO);
						LinkPath S0(gauge_field, T, L, n), ST(gauge_field, T, L, n);
						ST.move(0, WL_T);
						for (int r = 0; r < WL_r; ++r) {
							S0(i, true);
							ST(i, true);
						}

						complex curr_WL = complete_Wilson_loop(SO, S0.path, ST.path);
						co_pl_eq_co(&WL_avg, &curr_WL);
					}
				}
			}
		}
	}

	co_di_eq_re(&WL_avg, 3.0 * T / 2.0 * L * L * L);

	cout << scientific << setprecision(11) << showpos << WL_avg.re << " " << WL_avg.im << "\n";
}
