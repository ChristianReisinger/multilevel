#include <array>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_TWOLINK_OPERATOR_FUNCTIONS_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_TWOLINK_OPERATOR_FUNCTIONS_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

using twolink_operator_sig = void(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep);

inline int z_rel_dir(int dir, int z_dir) {
	return (dir + z_dir - 1) % 3 + 1;
}

void clover(double* result, const double* sub_gauge_field, int T, int L,
		const std::array<int, 4>& n, int mu, int nu, bool bar = false);

void C_x_I(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep, int mu, int nu, bool bar = false);

void I_x_C(double* result, const double* sub_gauge_field, int T, int L,
		int t, int x, int y, int z, int dir, int rsep, int mu, int nu, bool bar = false);

/**
 * function naming:
 * '_x_' computes the 'two-link' operator as direct product of the expression to its left with the one to its right
 * 'U' is a single link in temporal direction, at R=0 on the left and R=rsep on the right of '_x_'
 * 'Ex, Exbar, Bx, ..' are E/B fields ('bar' for HM factor) computed with clover, 'low'/'upp' indicate that only the lower/upper half clover is used
 * 'I' is the identity
 */

twolink_operator_sig U_x_U;
twolink_operator_sig  UU_x_UU;

twolink_operator_sig UU_x_UUEzlow;
twolink_operator_sig UU_x_EzuppUU;

twolink_operator_sig Ex_x_I;
twolink_operator_sig Exbar_x_I;
twolink_operator_sig Ey_x_I;
twolink_operator_sig Eybar_x_I;
twolink_operator_sig Ez_x_I;
twolink_operator_sig Ezbar_x_I;

twolink_operator_sig Bx_x_I;
twolink_operator_sig Bxbar_x_I;
twolink_operator_sig By_x_I;
twolink_operator_sig Bybar_x_I;
twolink_operator_sig Bz_x_I;
twolink_operator_sig Bzbar_x_I;

twolink_operator_sig I_x_Ex;
twolink_operator_sig I_x_Exbar;
twolink_operator_sig I_x_Ey;
twolink_operator_sig I_x_Eybar;
twolink_operator_sig I_x_Ez;
twolink_operator_sig I_x_Ezbar;

twolink_operator_sig I_x_Bx;
twolink_operator_sig I_x_Bxbar;
twolink_operator_sig I_x_By;
twolink_operator_sig I_x_Bybar;
twolink_operator_sig I_x_Bz;
twolink_operator_sig I_x_Bzbar;

twolink_operator_sig IU_x_IU;

}
}
}

#endif /* INCLUDE_TWOLINK_OPERATORS_HH_ */
