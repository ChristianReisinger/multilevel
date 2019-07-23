#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <array>

#include <linear_algebra.hh>
#include <geometry.hh>
#include <geometry2.hh>
#include <fields.hh>
#include <io.hh>

//compute [ [T(0)T(1)] [T(2)T'(3)] ] [ [T''(4)T(5)] [T(6)T(7)] ]
//with T' T'' such that E_z is inserted at t=4 as a clover

const int T = 20, L = 20;
const int SUN_N = 2;
const int SUN_elems = 2 * SUN_N * SUN_N;
const std::string config_prefix = "/scratch/mesonqcd/reisinger/gauge_fields/su2_b2.74_L20T20/confs/useable_confs/conf";

class LinkPath {

	double path[SUN_elems];

	LinkPath(double* gauge_field, const std::array<int, 4>& n) :
			gauge_field(gauge_field), n(n) {
		cm_eq_id(path);
	}

	LinkPath& operator()(int mu, bool pos) {
		double U[SUN_elems];
		cm_eq_cm(U, path);
		if (pos) {
			cm_eq_cm_ti_cm(path, U, gauge_field + ggi_n(n[0], n[1], n[2], n[3], 4, mu, T, L));
			n[mu]++;
		} else {
			n[mu]--;
			cm_eq_cm_ti_cm_dag(path, U, gauge_field + ggi_n(n[0], n[1], n[2], n[3], 4, mu, T, L));
		}
		return *this;
	}

	LinkPath& move(int mu, int dist) {
		n[mu] += dist;
		return *this;
	}

private:
	double* gauge_field;
	std::array<int, 4> n;
};

std::string config_filename(const std::vector<int>& tag) {
	std::ostringstream filename_oss;
	filename_oss << config_prefix;
	for(int i : tag)
		filename_oss << "." << i;
	return filename_oss.str();
}

/**
 * n : starting lattice site of the full operator
 * t_sub : initial time of the timeslice on which to compute the sublattice operator
 */
void sublattice_operator(const std::vector<int>& conf_tag, int t_sub, const std::array<int, 4>& n, int rsep, int dir,
		double* left_ret, double* right_ret) {

	//TODO t_sub / t_offset incorrect ?

	double* sub_gauge_field;
	Gauge_Field_Alloc(&sub_gauge_field, T, L);
	read_gauge_field(sub_gauge_field, config_filename(conf_tag).c_str(), T, L);


	LinkPath left(sub_gauge_field, n), right(sub_gauge_field, n);
	left.move(0, t_sub);
	cm_eq_cm(left_ret, left(0, true)(0, true).path);

	right.move(0, t_sub).move(dir, rsep);
	right(0, true)(0, true);

	if (t_sub == 2) {
		double clov[SUN_elems], clov_dag[SUN_elems];
		LinkPath p(sub_gauge_field, n);
		p.move(0, t_sub + 2).move(dir, rsep);
		cm_eq_cm(clov, p(dir, true)(0, false)(dir, false)(0, true).path);

		p = LinkPath(sub_gauge_field, n);
		p.move(0, t_sub + 2).move(dir, rsep);
		cm_pl_eq_cm(clov, p(0, false)(dir, false)(0, true)(dir, true).path);

		cm_eq_cm_dag(clov_dag, clov);
		cm_ti_eq_re(clov_dag, -1);
		cm_pl_eq_cm(clov, clov_dag);
		cm_ti_eq_re(clov, 0.5);
		cm_eq_cm_ti_cm(right_ret, right.path, clov);
	} else if (t_sub == 4) {
		double clov[SUN_elems], clov_dag[SUN_elems];
		LinkPath p(sub_gauge_field, n);
		cm_eq_cm(clov, p(0, true)(dir, true)(0, false)(dir, false).path);

		p = LinkPath(sub_gauge_field, n);
		cm_pl_eq_cm(clov, p(dir, false)(0, true)(dir, true)(0, false).path);

		cm_eq_cm_dag(clov_dag, clov);
		cm_ti_eq_re(clov_dag, -1);
		cm_pl_eq_cm(clov, clov_dag);
		cm_ti_eq_re(clov, 0.5);
		cm_eq_cm_ti_cm(right_ret, clov, right.path);
	} else {
		cm_eq_cm(right_ret, right.path);
	}
}

const int thickness_lv0 = 4;
const int updates_lv1 = 50;
const int thickness_lv1 = 2;
const int updates_lv2 = 20;

void multilevel(const std::vector<int>& conf_tag, int level, int t_offset, const std::array<int, 4>& n, int dir, int rsep,
		double* left_ret, double* right_ret) {
	if (level = 0) {
		cm_eq_id(left_ret);
		cm_eq_id(right_ret);

		for (int t_sub = 0; t_sub < T / thickness_lv0; t_sub += thickness_lv0) {
			double left_subavg[SUN_elems], right_subavg[SUN_elems], U[SUN_elems];
			multilevel(conf_tag, 1, t_sub, n, dir, rsep, left_subavg, right_subavg);
			cm_eq_cm_ti_cm(U, left_ret, left_subavg);
			cm_eq_cm(left_ret, U);
			cm_eq_cm_ti_cm(U, right_ret, right_subavg);
			cm_eq_cm(right_ret, U);
		}
	} else if (level = 1) {
		cm_eq_zero(left_ret);
		cm_eq_zero(right_ret);

		for (int i_lv1 = 0; i_lv1 < updates_lv1; ++i_lv1) {
			std::vector<int> curr_tag(conf_tag);
			curr_tag.push_back(i_lv1);

			double curr_left[SUN_elems], curr_right[SUN_elems];
			cm_eq_id(curr_left);
			cm_eq_id(curr_right);

			for (int t_sub = 0; t_sub < thickness_lv0 / thickness_lv1; t_sub += thickness_lv1) {
				double left_subavg[SUN_elems], right_subavg[SUN_elems], U[SUN_elems];
				multilevel(curr_tag, 2, t_offset + t_sub, n, dir, rsep, left_subavg, right_subavg);
				cm_eq_cm_ti_cm(U, curr_left, left_subavg);
				cm_eq_cm(curr_left, U);
				cm_eq_cm_ti_cm(U, curr_right, left_subavg);
				cm_eq_cm(curr_right, U);
			}

			cm_pl_eq_cm(left_ret, curr_left);
			cm_pl_eq_cm(right_ret, curr_right);
		}
		cm_ti_eq_re(left_ret, 1.0 / updates_lv1);
		cm_ti_eq_re(right_ret, 1.0 / updates_lv1);
	} else if (level = 2) {
		cm_eq_zero(left_ret);
		cm_eq_zero(right_ret);

		for (int i_lv2 = 0; i_lv2 < updates_lv2; ++i_lv2) {
			std::vector<int> curr_tag(conf_tag);
			curr_tag.push_back(i_lv2);

			double left_subop[SUN_elems], right_subop[SUN_elems];
			sublattice_operator(curr_tag, t_offset, n, rsep, dir, left_subop, right_subop);
			cm_pl_eq_cm(left_ret, left_subop);
			cm_pl_eq_cm(right_ret, right_subop);
		}
		cm_ti_eq_re(left_ret, 1.0 / updates_lv2);
		cm_ti_eq_re(right_ret, 1.0 / updates_lv2);
	}
}

int main(int argc, char **argv) {
	int WL_r = 5, WL_t = 8;

	if(argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <config_lv0_id>\n";
		return 0;
	}

	int config_lv0_id;
	std::stringstream config_id_ss;
	config_id_ss << argv[1];
	config_id_ss >> config_lv0_id;
	if(config_id_ss.fail()) {
		std::cerr << "Error: could not read <config_lv0_id>\n";
		return 0;
	}


	complex WL_avg;
	co_eq_zero(&WL_avg);

	for (int t = 0; t < T; t += 2) {
		for (int x = 0; x < L; ++x) {
			for (int y = 0; y < L; ++y) {
				for (int z = 0; z < L; ++z) {
					std::array<int, 4> n = { t, x, y, z };
					for (int i = 1; i < 4; ++i) {

						double left[SUN_elems], right[SUN_elems];
						multilevel( { config_lv0_id }, 0, 0, n, i, WL_r, left, right);
						LinkPath S0(n), St(n);
						St.move(0, WL_t);
						for (int r = 0; r < WL_r; ++r) {
							S0(i, true);
							St(i, true);
						}

						double U[SUN_elems], V[SUN_elems];
						cm_eq_cm_ti_cm(U, S0.path, right);
						cm_eq_cm_ti_cm_dag(V, U, St.path);
						cm_eq_cm_ti_cm_dag(U, V, left);
						complex WL;
						co_eq_tr_cm(&WL, U);
						WL.re /= SUN_N;
						WL.im /= SUN_N;
						co_pl_eq_co(&WL_avg, &WL);
					}
				}
			}
		}
	}

	std::cout << WL_avg.re / (3 * T / 2.0 * L * L * L);
}
