#include <sstream>
#include <chrono>
#include <ctime>
#include <iomanip>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {
namespace logger {

std::string timestamp() {
	std::time_t t = std::time(nullptr);
	std::ostringstream timestamp_oss;
	timestamp_oss << "[" << std::put_time(std::localtime(&t), "%T") << "] ";
	return timestamp_oss.str();
}

}
}
}
}
