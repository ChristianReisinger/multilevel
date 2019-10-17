#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <iomanip>
#include <cmath>
#include <chrono>

#include <fields.hh>
#include <io.hh>
#include <ranlux.hh>
#include <heatbath.hh>

#include <LevelDef.hh>

#include <MultilevelConfig.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

MultilevelConfig::MultilevelConfig(const std::string& filename_prefix, int top_level_id, int T, int L,
		const std::vector<LevelDef>& levels, double beta = 0, int seed = 0, bool save = false) :
		filename_prefix(filename_prefix), curr_tag( { top_level_id }), T(T), L(L),
				beta(beta), save(save), generate_configs(beta > 0 && seed > 0), time_spent_generating(0) {

	Gauge_Field_Alloc_silent(&top_level_conf, T, L);
	Gauge_Field_Alloc_silent(&config_buf, T, L);

	read_gauge_field(top_level_conf, config_filename().c_str(), T, L);

	for (const LevelDef& level : levels)
		m_levels.push_back(&level);

	if (generate_configs) {
		InitializeRand(seed);
		std::cerr << "Updating top level config ... ";
		for (int i_swp = 0; i_swp < levels.at(0).update_num(); ++i_swp)
			do_sweep(top_level_conf, T, L, beta);
		std::cerr << "ok\n";
	}
}

MultilevelConfig::~MultilevelConfig() {
	Gauge_Field_Free(&config_buf);
	Gauge_Field_Free(&top_level_conf);
}

int MultilevelConfig::config_num(int level) const {
	return m_levels.at(level)->config_num();
}

std::vector<int> MultilevelConfig::thickness(int level) const {
	return m_levels.at(level)->timeslice_sizes();
}

void MultilevelConfig::get(double*& gauge_field) const {
	gauge_field = config_buf;
}

std::string MultilevelConfig::config_filename() const {
	return filename_prefix + tag_to_string();
}

void MultilevelConfig::update(int level) {
	next_tag(level);
	if (level == 0) {
		Gauge_Field_Copy(config_buf, top_level_conf, T, L);
	} else if (generate_configs) {
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
			do_sweep(config_buf, T, L, beta, boundary_ts);

		time_spent_generating += std::chrono::steady_clock::now() - start_time;
		std::cerr << "ok\n";
	} else if (level == m_levels.size() - 1) {
		std::cerr << "Reading config '" << config_filename() << "' ... ";
		read_gauge_field(config_buf, config_filename().c_str(), T, L);
		std::cerr << "ok\n";
	}

	if (generate_configs && save)
		write_config();
}

int MultilevelConfig::milliseconds_spent_generating() const {
	return std::chrono::duration_cast<std::chrono::milliseconds>(time_spent_generating).count();
}

// private

void MultilevelConfig::next_tag(int level) {
	if (level == 0)
		curr_tag = std::vector<int> { curr_tag.at(0) };
	else if (level == curr_tag.size() - 1)
		++curr_tag.at(level);
	else if (level == curr_tag.size())
		curr_tag.push_back(1);
	else if (level < curr_tag.size() - 1) {
		curr_tag.erase(curr_tag.begin() + level + 1, curr_tag.end());
		++curr_tag.at(level);
	}
}

std::string MultilevelConfig::tag_to_string() const {
	std::ostringstream tag_oss;
	for (int i = 0; i < curr_tag.size(); ++i)
		tag_oss << "." << std::setfill('0') << std::setw(log10(m_levels.at(i)->config_num()) + 1) << curr_tag.at(i);
	return tag_oss.str();
}

void MultilevelConfig::write_config() const {
	std::string config_filename = filename_prefix + ".multilevel" + tag_to_string();
	std::ostringstream header_oss;
	header_oss << "generated during multilevel : " << beta << " " << T << " " << L;
	write_gauge_field(config_buf, config_filename.c_str(), T, L, header_oss.str().c_str());
}

}
}
}
