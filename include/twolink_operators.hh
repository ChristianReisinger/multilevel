#ifndef INCLUDE_TWOLINK_OPERATORS_HH_
#define INCLUDE_TWOLINK_OPERATORS_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

inline int z_rel_dir(int dir, int z_dir) {
	return (dir + z_dir - 1) % 3 + 1;
}

void clover(double* result, const double* sub_gauge_field, int T, int L,
		const std::array<int, 4>& n, int mu, int nu, bool bar = false);

void C_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep, int mu, int nu, bool bar = false);

void I_x_C(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep, int mu, int nu, bool bar = false);

/**
 * function naming:
 * '_x_' computes the 'two-link' operator as direct product of the expression to its left with the one to its right
 * 'U' is a single link in temporal direction, at R=0 on the left and R=rsep on the right of '_x_'
 * 'Ex, Exbar, Bx, ..' are E/B fields ('bar' for HM factor) computed with clover, 'low'/'upp' indicate that only the lower/upper half clover is used
 * 'I' is the identity
 */

void U_x_U(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void UU_x_UU(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void UU_x_UUEzlow(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void UU_x_EzuppUU(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void Ex_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Exbar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Ey_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Eybar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Ez_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Ezbar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Bx_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Bxbar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void By_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Bybar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Bz_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void Bzbar_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void I_x_Ex(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Exbar(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Ey(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Eybar(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Ez(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Ezbar(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Bx(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Bxbar(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_By(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Bybar(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Bz(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);
void I_x_Bzbar(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

void IU_x_IU(double* result, const double* sub_gauge_field, int T, int L,
		int& t, int x, int y, int z, int dir, int rsep);

}
}
}

#endif /* INCLUDE_TWOLINK_OPERATORS_HH_ */
