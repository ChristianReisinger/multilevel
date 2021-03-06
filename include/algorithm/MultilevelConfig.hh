#include <string>
#include <vector>
#include <chrono>
#include <memory>

#include <SUNGaugefield.hh>

#include <algorithm/LevelDef.hh>
#include <parameters/ConfigParameters.hh>

#ifndef INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_MULTILEVELCONFIG_HH_
#define INCLUDE_DE_UNI_FRANKFURT_ITP_REISINGER_MULTILEVEL_0819_MULTILEVELCONFIG_HH_

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace multilevel_0819 {

class MultilevelConfig {
public:

	MultilevelConfig(ConfigParameters config_params);

	~MultilevelConfig();
	MultilevelConfig(const MultilevelConfig&) = delete;
	MultilevelConfig(MultilevelConfig&&) = delete;
	MultilevelConfig& operator=(const MultilevelConfig&) = delete;
	MultilevelConfig& operator=(MultilevelConfig&&) = delete;

	std::string config_filepath() const;
	const double* get() const;
	int milliseconds_spent_updating() const;

	int get_T() const;
	int get_L() const;

private:
	void set_levels(const std::vector<LevelDef*>& levels);
	void update(size_t level);

	void next_tag(size_t level);
	std::string tag_to_string() const;
	void write_config() const;

	const std::string m_filestem;
	std::vector<const LevelDef*> m_levels { };

	const double m_beta;
	const int m_seed;
	const bool m_write;
	const bool m_generate;

	std::vector<int> m_tag;
	double* m_top_level_conf;

	std::unique_ptr<latticetools_0719::SUNGaugefield> m_SUN_gaugefield;

	std::chrono::steady_clock::duration m_time_spent_updating { 0 };

	friend class LevelAccess;
};

}
}
}

#endif /* INCLUDE_MULTILEVELCONFIG_HH_ */
