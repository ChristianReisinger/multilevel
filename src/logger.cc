#include <iostream>
#include <chrono>
#include <ctime>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {
namespace logger {

void print_timestamp() {
	const std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	std::cout << "[" << std::ctime(&time) << "] ";
}

}
}
}
}
