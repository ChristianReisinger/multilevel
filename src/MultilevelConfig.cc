/*
 * MultilevelConfig.cc
 *
 *  Created on: 28 Aug 2019
 *      Author: reisinger
 */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>

#include <fields.hh>
#include <io.hh>
#include <ranlux.hh>
#include <heatbath.hh>

#include <MultilevelConfig.hh>

MultilevelConfig::MultilevelConfig(const std::string& filename_prefix, int top_level_id, int T, int L,
		const std::vector<int>& level_thickness, const std::vector<int>& level_config_num,
		int beta, int seed, std::vector<int> level_updates, bool save) :
		filename_prefix(filename_prefix), curr_tag( { top_level_id }), T(T), L(L),
				level_thickness(level_thickness), level_config_num(level_config_num),
				beta(beta), level_updates(level_updates), save(save),
				generate_configs(beta > 0 && seed > 0 && !level_updates.empty()), time_spent_generating(0) {

	Gauge_Field_Alloc_silent(&top_level_conf, T, L);
	Gauge_Field_Alloc_silent(&config_buf, T, L);

	read_gauge_field(top_level_conf, config_filename().c_str(), T, L);
	if (generate_configs) {
		InitializeRand(seed);
		for (int i_swp = 0; i_swp < level_updates[0]; ++i_swp)
			do_sweep(top_level_conf, T, L, beta);
	}
}

MultilevelConfig::~MultilevelConfig() {
	Gauge_Field_Free(&config_buf);
	Gauge_Field_Free(&top_level_conf);
}

int MultilevelConfig::thickness(int level) const {
	return level_thickness[level];
}

int MultilevelConfig::config_num(int level) const {
	return level_config_num[level];
}

void MultilevelConfig::get(double*& gauge_field) const {
	gauge_field = config_buf;
}

std::string MultilevelConfig::config_filename() const {
	std::ostringstream filename_oss;
	filename_oss << filename_prefix << tag_to_string();
	return filename_oss.str();
}

void MultilevelConfig::update(int level) {
	next_tag(level);
	if (level == 0) {
		Gauge_Field_Copy(config_buf, top_level_conf, T, L);
	} else if (generate_configs) {
		std::cerr << "Generating config '" << config_filename() << "' ... ";
		auto start_time = std::chrono::steady_clock::now();

		std::vector<int> boundary_ts;
		for (int t = 0; t < T; t += level_thickness[level])
			boundary_ts.push_back(t);
		for (int i_swp = 0; i_swp < level_updates[level]; ++i_swp)
			do_sweep(config_buf, T, L, beta, boundary_ts);

		time_spent_generating += std::chrono::steady_clock::now() - start_time;
		std::cerr << "ok\n";
	} else if (level == level_config_num.size() - 1) {
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
		curr_tag = std::vector<int> { curr_tag[0] };
	else if (level == curr_tag.size() - 1)
		++curr_tag[level];
	else if (level == curr_tag.size())
		curr_tag.push_back(1);
	else if (level < curr_tag.size() - 1) {
		curr_tag.erase(curr_tag.begin() + level + 1, curr_tag.end());
		++curr_tag[level];
	}
}

std::string MultilevelConfig::tag_to_string() const {
	std::ostringstream tag_oss;
	for (int i = 0; i < curr_tag.size(); ++i)
		tag_oss << "." << std::setfill('0') << std::setw(log10(level_config_num.at(i)) + 1) << curr_tag.at(i);
	return tag_oss.str();
}

void MultilevelConfig::write_config() const {
	std::ostringstream config_filename_oss;
	config_filename_oss << filename_prefix << ".multilevel" << tag_to_string();
	std::ostringstream header_oss;
	header_oss << "generated during multilevel : " << beta << " " << T << " " << L;
	write_gauge_field(config_buf, config_filename_oss.str().c_str(), T, L, header_oss.str().c_str());
}
