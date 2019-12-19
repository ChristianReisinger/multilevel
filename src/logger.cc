#include <iostream>
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

}
}
}
}
