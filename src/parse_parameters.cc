#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

#include <global_defs.hh>
#include <linear_algebra.hh>
#include <helper_functions.hh>

#include <TwolinkComputer.hh>
#include <LevelDef.hh>
#include <TwolinkOperator.hh>

#include <parse_parameters.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

namespace parse_parameters {

namespace {

static const std::regex format(""
		"(thickness \\d+(,\\d+)*\n+"
		"((?!thickness)\\S+?:[.x]+:(( |\t)+\\S+)+\n+)+)+"
		"thickness T\n+"
		"((?!thickness)\\S+?:[.x]+:.*:\\d+:(( |\t)+\\S+)+\n+)+", std::regex::nosubs);

static const std::regex level_format(""
		"thickness (\\d+(?:,\\d+)*|T)\n+"
		"((?:(?!thickness)\\S+?:[.x]+(?::.*:\\d+)?:(?:(?: |\t)+\\S+)+\n+)+)");

static const std::string name_tsldef_format = "(\\S+?):([.x]+):";
static const std::string descr_T_format = "(.*):(\\d+):"; //TODO T is not needed anymore .. computed automatically in TwolinkOperator class
static const std::string factors_regex = "((?:(?: |\t)+\\S+)+)\n";
static const std::regex operator_format(name_tsldef_format + factors_regex);
static const std::regex operator_format_top(name_tsldef_format + descr_T_format + factors_regex);

template<typename T>
std::vector<const OperatorFactor*> parse_operator_factors(const std::vector<T>& available_operators, const std::string& factor_str) {
	std::vector<const OperatorFactor*> factors;
	static const std::regex factor_format("(\\S+)");
	for (std::sregex_iterator factor_it(factor_str.begin(), factor_str.end(), factor_format);
			factor_it != std::sregex_iterator(); ++factor_it) {

		const std::string requested_factor_name = factor_it->str(1);
		const auto& requested_factor_it = std::find_if(available_operators.begin(), available_operators.end(), [&](const T& f) {
			return f.name() == requested_factor_name;
		});

		if (requested_factor_it != available_operators.end()) {
			const OperatorFactor* f = &*requested_factor_it;
			factors.push_back(f);
		} else
			throw std::invalid_argument("Twolink operator '" + requested_factor_name + "' does not exist");
	}
	return factors;
}

std::vector<bool> parse_timeslice_defined(const std::string& timeslice_def_str) {
	std::vector<bool> timeslice_defined;
	std::istringstream is_def_iss(timeslice_def_str);
	char is_def_char;
	while (is_def_iss >> is_def_char) {
		if (is_def_char == '.')
			timeslice_defined.push_back(false);
		else if (is_def_char == 'x')
			timeslice_defined.push_back(true);
	}
	return timeslice_defined;
}

std::pair<bool, std::regex> get_operator_format(const std::string& operators_str) {
	if (std::regex_search(operators_str, operator_format_top))
		return {true, operator_format_top};
	else
		return {false, operator_format};
}

template<typename T>
void parse_operators(LevelDef& level, const std::vector<T>& available_operators, const std::string& operators_str) {

	auto istop_fmt = get_operator_format(operators_str);
	for (std::sregex_iterator operator_it(operators_str.begin(), operators_str.end(), istop_fmt.second);
			operator_it != std::sregex_iterator(); ++operator_it) {

		const std::string operator_name = operator_it->str(1);

		TwolinkOperator curr_op(operator_name,
				parse_timeslice_defined(operator_it->str(2)),
				parse_operator_factors(available_operators, operator_it->str(istop_fmt.first ? 5 : 3)));

		if (istop_fmt.first)
			curr_op.descr(operator_it->str(3));

		if (istop_fmt.first || std::count_if(level.operators().begin(), level.operators().end(), [&](const TwolinkOperator& op) {
			return op.name() == operator_name;
		}) == 0)
			level.add_operator(curr_op);
		else
			throw std::runtime_error("duplicate definition of '" + operator_name + "'");
	}
}

bool verify_timeslice_sizes(const std::vector<LevelDef>& levels) {
	if (levels.empty())
		return false;
	for (auto level_it = ++levels.begin(); level_it != levels.end(); ++level_it) {

		const auto& tsl_sizes = level_it->timeslice_sizes();
		auto size_it = tsl_sizes.begin();
		for (const int prev_level_size : std::prev(level_it)->timeslice_sizes()) {
			int curr_total_size = 0;

			while (size_it != tsl_sizes.end()) {
				curr_total_size += *size_it;
				++size_it;
				if (curr_total_size == prev_level_size)
					break;
				else if (curr_total_size > prev_level_size)
					return false;
			}
			if (curr_total_size != prev_level_size)
				return false;
		}
		if (size_it != tsl_sizes.end())
			return false;
	}
	return true;
}

}

std::vector<LevelDef> levels(const std::vector<TwolinkComputer>& twolink_computers, const std::string& compstr) {
	std::vector<LevelDef> levels;
	if (!std::regex_match(compstr, format))
		throw std::runtime_error("invalid composition format");

	const std::sregex_iterator level_begin(compstr.begin(), compstr.end(), level_format);
	if (std::distance(level_begin, std::sregex_iterator()) < 2)
		throw std::invalid_argument("less than 2 levels");

	for (auto level_it = level_begin; level_it != std::sregex_iterator(); ++level_it) {

		std::vector<int> timeslice_sizes;
		const std::string size_str = level_it->str(1);
		if (size_str == "T")
			timeslice_sizes = levels[0].timeslice_sizes();
		else
			timeslice_sizes = tools::helper::parse_unsigned_int_list(size_str.c_str());
		levels.insert(levels.begin(), LevelDef(timeslice_sizes));

		if (level_it == level_begin)
			parse_operators(levels[0], twolink_computers, level_it->str(2));
		else
			parse_operators(levels[0], levels[1].operators(), level_it->str(2));
	}
	if (!verify_timeslice_sizes(levels))
		throw std::invalid_argument("invalid timeslice sizes");

	return levels;
}

}
}
}
}
