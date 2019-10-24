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
#include <fields.hh>
#include <io.hh>
#include <ranlux.hh>
#include <heatbath.hh>

#include <LevelDef.hh>

#include <MultilevelConfig.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

MultilevelConfig::MultilevelConfig(const std::string& filestem, int top_level_id, int T, int L,
		double beta, int seed, bool write) :
		m_filestem(filestem), m_tag( { top_level_id }), T(T), L(L),
				m_beta(beta), m_seed(seed), m_write(write), m_generate(beta > 0 && seed > 0) {

	if (filestem.empty() || top_level_id < 0)
		throw std::invalid_argument("invalid MultilevelConfig");

	Gauge_Field_Alloc_silent(&m_top_level_conf, T, L);
	Gauge_Field_Alloc_silent(&m_config_buf, T, L);
}

MultilevelConfig::~MultilevelConfig() {
	Gauge_Field_Free(&m_config_buf);
	Gauge_Field_Free(&m_top_level_conf);
}

void MultilevelConfig::get(double*& gauge_field) const {
	gauge_field = m_config_buf;
}

std::string MultilevelConfig::config_filename() const {
	return m_filestem + tag_to_string();
}

int MultilevelConfig::milliseconds_spent_generating() const {
	return std::chrono::duration_cast<std::chrono::milliseconds>(m_time_spent_generating).count();
}

// private

void MultilevelConfig::set_levels(std::vector<LevelDef*> levels) {
	for (const auto* level : levels)
		m_levels.push_back(level);

	read_gauge_field(m_top_level_conf, config_filename().c_str(), T, L);

	if (m_generate) {
		InitializeRand(m_seed);
		std::cerr << "Updating top level config ... ";
		for (int i_swp = 0; i_swp < m_levels[0]->update_num(); ++i_swp)
			do_sweep(m_top_level_conf, T, L, m_beta);
		std::cerr << "ok\n";
	}
}

void MultilevelConfig::update(int level) {
	next_tag(level);
	if (level == 0) {
		Gauge_Field_Copy(m_config_buf, m_top_level_conf, T, L);
	} else if (m_generate) {
		std::cerr << "Generating config '" << config_filename() << "' ... ";
		auto start_time = std::chrono::steady_clock::now();

		std::set<int> boundary_ts;
		int boundary_t = 0;
		while (boundary_t < T) {
			for (int timeslice_size : m_levels.at(level)->timeslice_sizes()) {
				boundary_ts.insert(boundary_t);
				boundary_t += timeslice_size;
			}
		}
		for (int i_swp = 0; i_swp < m_levels.at(level)->update_num(); ++i_swp)
			do_sweep(m_config_buf, T, L, m_beta, boundary_ts);

		m_time_spent_generating += std::chrono::steady_clock::now() - start_time;
		std::cerr << "ok\n";
	} else if (level == m_levels.size() - 1) {
		std::cerr << "Reading config '" << config_filename() << "' ... ";
		read_gauge_field(m_config_buf, config_filename().c_str(), T, L);
		std::cerr << "ok\n";
	}

	if (m_generate && m_write)
		write_config();
}

void MultilevelConfig::next_tag(int level) {
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
	for (int level = 0; level < m_tag.size(); ++level)
		tag_oss << "." << std::setfill('0') << std::setw(log10(m_levels.at(level)->config_num()) + 1) << m_tag.at(level);
	return tag_oss.str();
}

void MultilevelConfig::write_config() const {
	std::string config_filename = m_filestem + ".multilevel" + tag_to_string();
	std::ostringstream header_oss;
	header_oss << "generated during multilevel : " << m_beta << " " << T << " " << L;
	write_gauge_field(m_config_buf, config_filename.c_str(), T, L, header_oss.str().c_str());
}

}
}
}
