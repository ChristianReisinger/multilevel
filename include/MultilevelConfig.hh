/*
 * MultilevelConfig.hh
 *
 *  Created on: 28 Aug 2019
 *      Author: reisinger
 */

#include <string>
#include <vector>
#include <chrono>

#ifndef INCLUDE_MULTILEVELCONFIG_HH_
#define INCLUDE_MULTILEVELCONFIG_HH_

class MultilevelConfig {
public:

	MultilevelConfig(const std::string& filename_prefix, int top_level_id, int T, int L,
			const std::vector<int>& level_thickness, const std::vector<int>& level_config_num,
			double beta = 0, int seed = 0, std::vector<int> level_updates = { }, bool save = false);

	~MultilevelConfig();
	MultilevelConfig(const MultilevelConfig&) = delete;
	MultilevelConfig(MultilevelConfig&&) = delete;
	MultilevelConfig& operator=(const MultilevelConfig&) = delete;
	MultilevelConfig& operator=(MultilevelConfig&&) = delete;

	void update(int level);

	const int T, L;
	std::string config_filename() const;
	int thickness(int level) const;
	int config_num(int level) const;
	void get(double*& gauge_field) const;
	int milliseconds_spent_generating() const;

private:

	void next_tag(int level);
	std::string tag_to_string() const;
	void write_config() const;

	const std::string filename_prefix;
	const std::vector<int> level_thickness;
	const std::vector<int> level_config_num;

	const bool generate_configs;
	const double beta;
	const bool save;
	const std::vector<int> level_updates;

	std::vector<int> curr_tag;
	double* top_level_conf;
	double* config_buf;

	std::chrono::steady_clock::duration time_spent_generating;
};

#endif /* INCLUDE_MULTILEVELCONFIG_HH_ */
