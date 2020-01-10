#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <stdexcept>

#include <global_defs.hh>
#include <helper_functions.hh>
#if __SUN_N__ == 2
#include <MCSU2Gaugefield.hh>
#elif __SUN_N__ == 3
#include <CL2QCDGaugefield.hh>
#else
#	error INVALID NC
#endif

#include <LevelDef.hh>
#include <logger.hh>

#include <MultilevelConfig.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

MultilevelConfig::MultilevelConfig(const std::string& filestem, int top_level_id, int T, int L,
		double beta, int seed, int overrelax_steps, bool write) :
		m_filestem(filestem),
				m_beta(beta), m_seed(seed), m_write(write), m_generate(beta > 0 && seed > 1),
				m_tag( { top_level_id }) {

	if (filestem.empty() || top_level_id < 0)
		throw std::invalid_argument("invalid MultilevelConfig");

#if __SUN_N__ == 2
	m_SUN_gaugefield = tools::helper::make_unique<latticetools_0719::MCSU2Gaugefield>(T, L, seed, beta, config_filepath());
#elif __SUN_N__ == 3
	m_SUN_gaugefield = tools::helper::make_unique<latticetools_0719::CL2QCDGaugefield>(T, L, seed, beta, overrelax_steps, config_filepath());
#endif

	latticetools_0719::Gauge_Field_Alloc(m_top_level_conf, T, L);
}

MultilevelConfig::~MultilevelConfig() {
	latticetools_0719::Gauge_Field_Free(m_top_level_conf);
}

const double* MultilevelConfig::get() const {
	return m_SUN_gaugefield->get();
}

std::string MultilevelConfig::config_filepath() const {
	return m_filestem + tag_to_string();
}

int MultilevelConfig::milliseconds_spent_updating() const {
	return std::chrono::duration_cast<std::chrono::milliseconds>(m_time_spent_updating).count();
}

int MultilevelConfig::get_T() const {
	return m_SUN_gaugefield->get_T();
}

int MultilevelConfig::get_L() const {
	return m_SUN_gaugefield->get_L();
}

// private

void MultilevelConfig::set_levels(const std::vector<LevelDef*>& levels) {
	m_levels = std::vector<const LevelDef*>(levels.begin(), levels.end());

	if (m_generate) {
		std::cout << logger::timestamp() << "Updating top level config ... " << std::endl;
		for (int i_swp = 0; i_swp < m_levels[0]->update_num(); ++i_swp)
			m_SUN_gaugefield->do_sweep();
		std::cout << logger::timestamp() << "done\n";
	}
	latticetools_0719::Gauge_Field_Copy(m_top_level_conf, m_SUN_gaugefield->get(), get_T(), get_L());
}

void MultilevelConfig::update(size_t level) {
	auto start_time = std::chrono::steady_clock::now();

	next_tag(level);
	if (level == 0) {
		m_SUN_gaugefield->set(m_top_level_conf);
	} else if (m_generate) {
		std::cout << logger::timestamp() << "Generating config '" << config_filepath() << "' ... " << std::endl;

		std::set<int> fixed_timeslices;
		int boundary_t = 0;
		while (boundary_t < get_T()) {
			for (int timeslice_size : m_levels.at(level)->timeslice_sizes()) {
				fixed_timeslices.insert(boundary_t);
				boundary_t += timeslice_size;
			}
		}
		for (int i_swp = 0; i_swp < m_levels.at(level)->update_num(); ++i_swp)
			m_SUN_gaugefield->do_sweep(fixed_timeslices);

		std::cout << logger::timestamp() << "done\n";
	} else if (level == m_levels.size() - 1) {
		std::cout << logger::timestamp() << "Reading config '" << config_filepath() << "' ... " << std::endl;
		m_SUN_gaugefield->read(config_filepath());
		std::cout << logger::timestamp() << "done\n";
	}

	if (m_generate && m_write)
		write_config();

	m_time_spent_updating += std::chrono::steady_clock::now() - start_time;
}

void MultilevelConfig::next_tag(size_t level) {
	if (level == 0)
		m_tag = std::vector<int> { m_tag.at(0) };
	else if (level == m_tag.size() - 1)
		++m_tag.at(level);
	else if (level == m_tag.size())
		m_tag.push_back(1);
	else if (level < m_tag.size() - 1) {
		m_tag.erase(m_tag.begin() + level + 1, m_tag.end());
		++m_tag.at(level);
	}
}

std::string MultilevelConfig::tag_to_string() const {
	std::ostringstream tag_oss;
	for (size_t level = 0; level < m_tag.size(); ++level)
		tag_oss << "." << m_tag.at(level);
	return tag_oss.str();
}

void MultilevelConfig::write_config() const {
	std::string config_filename = m_filestem + ".multilevel" + tag_to_string();
	m_SUN_gaugefield->write(config_filename);
}

}
}
}
