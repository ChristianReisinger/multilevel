#include <array>

#include <LinkPath.hh>
#include <linear_algebra.hh>
#include <sublattice_algebra.hh>
#include <twolink_operator_functions.hh>
#include <iostream>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

using latticetools_0719::SUN_elems;
using latticetools_0719::LinkPath;

void clover(double* result, const double* sub_gauge_field, int T, int L,
		const std::array<int, 4>& n, int mu, int nu, bool bar) {

	LinkPath plaq(sub_gauge_field, T, L, n);

	plaq(mu, true)(nu, true)(mu, false)(nu, false);
	cm_eq_cm(result, plaq.path);

	plaq.reset(n)(nu, false)(mu, true)(nu, true)(mu, false);
	cm_pl_eq_cm(result, plaq.path);

	plaq.reset(n)(mu, false)(nu, false)(mu, true)(nu, true);
	cm_pl_eq_cm(result, plaq.path);

	plaq.reset(n)(nu, true)(mu, false)(nu, false)(mu, true);
	cm_pl_eq_cm(result, plaq.path);

	double U[SUN_elems];
	cm_eq_cm_dag(U, result);
	if (!bar)
		cm_ti_eq_re(U, -1.0);
	cm_pl_eq_cm(result, U);
	cm_ti_eq_re(result, 0.125);
}

void C_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep, int mu, int nu, bool bar) {
	double clov[SUN_elems];
	clover(clov, sub_gauge_field, T, L, { t, x, y, z }, mu, nu, bar);

	double I[SUN_elems];
	cm_eq_id(I);

	so_eq_cm_x_cm(result, clov, I);
}

void I_x_C(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep, int mu, int nu, bool bar) {
	double clov[SUN_elems];
	std::array<int, 4> n_clov = { t, x, y, z };
	n_clov[dir] += rsep;
	clover(clov, sub_gauge_field, T, L, n_clov, mu, nu, bar);

	double I[SUN_elems];
	cm_eq_id(I);

	so_eq_cm_x_cm(result, I, clov);
}

void U_x_U(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true);

	so_eq_cm_x_cm(result, T0.path, TR.path);
}

void UU_x_UU(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true)(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true)(0, true);

	so_eq_cm_x_cm(result, T0.path, TR.path);
}

void UU_x_UUEzlow(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true)(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true)(0, true);

	double clov[SUN_elems];
	LinkPath plaq(sub_gauge_field, T, L, { t + 2, x, y, z });

	plaq.move(dir, rsep)(dir, true)(0, false)(dir, false)(0, true);
	cm_eq_cm(clov, plaq.path);

	plaq.reset( { t + 2, x, y, z }).move(dir, rsep)(0, false)(dir, false)(0, true)(dir, true);
	cm_pl_eq_cm(clov, plaq.path);

	double U[SUN_elems];
	cm_eq_cm_dag(U, clov);
	cm_ti_eq_re(U, -1.0);
	cm_pl_eq_cm(clov, U);
	cm_ti_eq_re(clov, 0.5);

	cm_eq_cm_ti_cm(U, TR.path, clov);

	so_eq_cm_x_cm(result, T0.path, U);
}

void UU_x_EzuppUU(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t, x, y, z });
	T0(0, true)(0, true);

	LinkPath TR(sub_gauge_field, T, L, { t, x, y, z });
	TR.move(dir, rsep)(0, true)(0, true);

	double clov[SUN_elems];
	LinkPath plaq(sub_gauge_field, T, L, { t, x, y, z });

	plaq.move(dir, rsep)(0, true)(dir, true)(0, false)(dir, false);
	cm_eq_cm(clov, plaq.path);

	plaq.reset( { t, x, y, z }).move(dir, rsep)(dir, false)(0, true)(dir, true)(0, false);
	cm_pl_eq_cm(clov, plaq.path);

	double U[SUN_elems];
	cm_eq_cm_dag(U, clov);
	cm_ti_eq_re(U, -1.0);
	cm_pl_eq_cm(clov, U);
	cm_ti_eq_re(clov, 0.5);

	cm_eq_cm_ti_cm(U, clov, TR.path);

	so_eq_cm_x_cm(result, T0.path, U);
}

void I_x_Ex(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, z_rel_dir(1, dir));
}
void I_x_Exbar(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, z_rel_dir(1, dir), true);
}
void I_x_Ey(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, z_rel_dir(2, dir));
}
void I_x_Eybar(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, z_rel_dir(2, dir), true);
}
void I_x_Ez(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, dir);
}
void I_x_Ezbar(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, dir, true);
}
void I_x_Bx(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(2, dir), z_rel_dir(3, dir));
}
void I_x_Bxbar(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(2, dir), z_rel_dir(3, dir), true);
}
void I_x_By(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(3, dir), z_rel_dir(1, dir));
}
void I_x_Bybar(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(3, dir), z_rel_dir(1, dir), true);
}
void I_x_Bz(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(1, dir), z_rel_dir(2, dir));
}
void I_x_Bzbar(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	I_x_C(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(1, dir), z_rel_dir(2, dir), true);
}

void Ex_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, z_rel_dir(1, dir));
}
void Exbar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, z_rel_dir(1, dir), true);
}
void Ey_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, z_rel_dir(2, dir));
}
void Eybar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, z_rel_dir(2, dir), true);
}
void Ez_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, dir);
}
void Ezbar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, 0, dir, true);
}
void Bx_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(2, dir), z_rel_dir(3, dir));
}
void Bxbar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(2, dir), z_rel_dir(3, dir), true);
}
void By_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(3, dir), z_rel_dir(1, dir));
}
void Bybar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(3, dir), z_rel_dir(1, dir), true);
}
void Bz_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(1, dir), z_rel_dir(2, dir));
}
void Bzbar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	C_x_I(result, sub_gauge_field, T, L, t, x, y, z, dir, rsep, z_rel_dir(1, dir), z_rel_dir(2, dir), true);
}

void IU_x_IU(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep) {
	LinkPath T0(sub_gauge_field, T, L, { t + 1, x, y, z });
	T0(0, true);
	LinkPath TR(sub_gauge_field, T, L, { t + 1, x, y, z });
	TR.move(dir, rsep)(0, true);

	so_eq_cm_x_cm(result, T0.path, TR.path);
}

}
}
}
