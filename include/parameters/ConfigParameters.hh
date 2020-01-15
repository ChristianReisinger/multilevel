#include <string>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_CONFIGPARAMETERS_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_CONFIGPARAMETERS_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

struct ConfigParameters {
	int T = 4, L = 4, config_lv0_id = 1, seed = 1, overrelax_steps = 0;
	double beta = 0.0;
	bool write = false;
	std::string filestem = "conf";

	inline int spatial_volume() {
		return L * L * L;
	}
};

}
}
}

#endif
