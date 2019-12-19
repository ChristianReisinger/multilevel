#include <iostream>
#include <chrono>
#include <ctime>
#include <iomanip>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {
namespace logger {

void print_timestamp() {
	std::time_t t = std::time(nullptr);
	std::cout << "[" << std::put_time(std::localtime(&t), "%T") << "] ";
}

unsigned long get_ms_since_and_reset(std::chrono::steady_clock::time_point& time) {
	auto now = std::chrono::steady_clock::now();
	auto ms_since = std::chrono::duration_cast<std::chrono::milliseconds>(now - time).count();
	time = now;
	return ms_since;
}

}
}
}
}
