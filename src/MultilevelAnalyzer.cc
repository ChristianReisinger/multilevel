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

#include <global_defs.hh>

#include <sublattice_algebra.hh>
#include <T_field.hh>
#include <LevelDef.hh>
#include <LevelAccess.hh>
#include <TwolinkOperatorWriter.hh>
#include <MultilevelConfig.hh>

#include <MultilevelAnalyzer.hh>

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

int MultilevelAnalyzer::milliseconds_spent_computing() const {
	return std::chrono::duration_cast<std::chrono::milliseconds>(time_spent_computing_operators).count();
}

//private

void MultilevelAnalyzer::compute_sublattice_fields(const size_t level) {
	const bool is_lowest = (level == m_levels.size() - 1);
	const int config_num = (level == 0 ? 1 : m_levels[level]->config_num());

	for (int conf = 1; conf <= config_num; ++conf) {
		LevelAccess::update(*m_config, level);

		const double* lowest_level_gauge_field = nullptr;
		if (is_lowest)
			lowest_level_gauge_field = m_config->get();
		else {
			std::cerr << "Allocating sublattice fields on config '" << m_config->config_filepath() << "' ... ";
			LevelAccess::alloc_operators(*m_levels[level + 1], m_WL_Rs, m_config->get_T(), m_config->get_L());
			std::cerr << "ok\n";
			compute_sublattice_fields(level + 1);
		}

		auto start_time = std::chrono::steady_clock::now();
		std::cerr << "Computing observables on config '" << m_config->config_filepath() << "' ... " << std::flush;
		for (auto& op : LevelAccess::operators(*m_levels[level])) {
			for (const int WL_R : m_WL_Rs) {
				try {
					for (const int t : op.defined_ts(WL_R)) {
#pragma omp parallel for collapse(3)
						for (int x = 0; x < m_config->get_L(); ++x) {
							for (int y = 0; y < m_config->get_L(); ++y) {
								for (int z = 0; z < m_config->get_L(); ++z) {
									for (int i = 1; i < 4; ++i) {
										double curr_operator[SO_elems];
										so_eq_id(curr_operator);
										int curr_t = t;

										for (const auto& factor : op.factors) {
											double T_factor[SO_elems];
											factor->at(T_factor, curr_t, x, y, z, i, WL_R,
													lowest_level_gauge_field, m_config->get_T(), m_config->get_L());
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
		std::cerr << "o.k.\n";
		time_spent_computing_operators += std::chrono::steady_clock::now() - start_time;

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
