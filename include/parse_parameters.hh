#include <string>
#include <vector>

#include <TwolinkComputer.hh>
#include <LevelDef.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_PARSE_PARAMETERS_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_PARSE_PARAMETERS_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

namespace parse_parameters {
	std::vector<LevelDef> levels(const std::vector<TwolinkComputer>& twolink_computers, const std::string& compstr);
}

}
}
}

#endif
