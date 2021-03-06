#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>
#include <chrono>
#include <stdexcept>
#include <algorithm>
#include <omp.h>

#include <Stopwatch.hh>
#include <global_defs.hh>

#include <algorithm/LevelAccess.hh>
#include <algorithm/LevelDef.hh>
#include <algorithm/MultilevelAnalyzer.hh>
#include <algorithm/MultilevelConfig.hh>
#include <algorithm/T_field.hh>
#include <algorithm/TwolinkOperatorWriter.hh>
#include <log/logger.hh>
#include <physics/sublattice_algebra.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

//public

MultilevelAnalyzer::MultilevelAnalyzer(std::vector<LevelDef>& levels, MultilevelConfig& multilevel_config, std::set<int> WL_Rs) :
		m_config(&multilevel_config), m_WL_Rs(WL_Rs) {

	auto is_negative = [](int i) {
		return i < 0;
	};
	if (levels.size() < 1 || WL_Rs.empty() || std::count_if(WL_Rs.begin(), WL_Rs.end(), is_negative))
		throw std::invalid_argument("invalid MultilevelAnalyzer");

	for (auto& level : levels)
		m_levels.push_back(&level);

	LevelAccess::set_levels(multilevel_config, m_levels);
}

void MultilevelAnalyzer::compute_T_fields() {
	LevelAccess::alloc_operators(*m_levels[0], m_WL_Rs, m_config->get_T(), m_config->get_L());
	compute_sublattice_fields(0);
	LevelAccess::update(*m_config, 0);
}

std::chrono::milliseconds::rep MultilevelAnalyzer::milliseconds_spent_computing() const {
	return time_spent_computing_operators.count();
}

//private

void MultilevelAnalyzer::compute_sublattice_fields(const size_t level) {
	const bool is_lowest = (level == m_levels.size() - 1);
	const int config_num = (level == 0 ? 1 : m_levels[level]->config_num());
	const int config_T = m_config->get_T();
	const int config_L = m_config->get_L();

	for (int conf = 1; conf <= config_num; ++conf) {
		LevelAccess::update(*m_config, level);
		const std::string curr_config_filepath = m_config->config_filepath();

		const double* lowest_level_gauge_field = nullptr;
		if (is_lowest)
			lowest_level_gauge_field = m_config->get();
		else {
			std::cout << logger::timestamp() << "Allocating sublattice fields on config '" << curr_config_filepath << "' ... ";
			LevelAccess::alloc_operators(*m_levels[level + 1], m_WL_Rs, config_T, config_L);
			std::cout << "done\n";
			compute_sublattice_fields(level + 1);
		}

		tools::Stopwatch compute_watch;
		std::cout << logger::timestamp() << "Computing observables on config '" << curr_config_filepath << "' ... " << std::endl;
		for (auto& op : LevelAccess::operators(*m_levels[level])) {
			for (const int WL_R : m_WL_Rs) {
				try {
					for (const int t : op.defined_ts(WL_R)) {
#pragma omp parallel for collapse(3)
						for (int x = 0; x < config_L; ++x) {
							for (int y = 0; y < config_L; ++y) {
								for (int z = 0; z < config_L; ++z) {
									for (int i = 1; i < 4; ++i) {
										double curr_operator[SO_elems];
										so_eq_id(curr_operator);
										int curr_t = t;

										for (const auto& factor : op.factors) {
											double T_factor[SO_elems];
											factor->at(T_factor, curr_t, x, y, z, i, WL_R,
													lowest_level_gauge_field, config_T, config_L);
											curr_t += factor->t_extent();

											double T_temp[SO_elems];
											so_eq_so_ti_so(T_temp, curr_operator, T_factor);
											so_eq_so(curr_operator, T_temp);
										}

										so_pl_eq_so(TwolinkOperatorWriter::field(op, WL_R).T_at(t, x, y, z, i), curr_operator);
									}
								}
							}
						}
					}
				} catch (std::out_of_range& e) {
					throw std::runtime_error("invalid definition of operator '" + op.name() + "'");
				}
			}
		}
		std::cout << logger::timestamp() << "done\n";
		time_spent_computing_operators += compute_watch.check();

		if (!is_lowest)
			LevelAccess::free_operators(*m_levels[level + 1]);
	}
	for (auto& op : LevelAccess::operators(*m_levels[level]))
		for (const int WL_R : m_WL_Rs)
			TwolinkOperatorWriter::field(op, WL_R) /= config_num;
}

}
}
}
