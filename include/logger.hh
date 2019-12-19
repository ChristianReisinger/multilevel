#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_LOGGER_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_LOGGER_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {
namespace logger {

void print_timestamp();
unsigned long get_ms_since_and_reset(std::chrono::steady_clock::time_point& time);

}
}
}
}

#endif
