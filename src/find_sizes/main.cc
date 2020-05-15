#include <iostream>
#include <set>
#include <map>
#include <numeric>

#include <io_tools.hh>
#include <helper_functions.hh>

#include <find_sizes/main.hh>

namespace de_uni_frankfurt_itp::reisinger::multilevel_0819::find_sizes {

namespace {

// ------------- add more modes here ----------------
auto get_boundary_forbidden_time_coords(t_extent operator_t_extent, observable_mode mode) {

	std::set<std::set<t_coord> > boundary_forbidden_t_coords;

	switch (mode) {
		case center:
			std::set<t_extent> center_coords = { operator_t_extent / 2 };
			if (operator_t_extent % 2 != 0)
				center_coords.insert(operator_t_extent / 2 + 1);
			boundary_forbidden_t_coords.insert(center_coords);
			break;
	}
	return boundary_forbidden_t_coords;
}

int next_timeslice_boundary(int t, const timeslice_sizes& timeslice_sizes_pattern) {
	const int tsl_num = timeslice_sizes_pattern.size();
	int curr_t = 0;

	for (int tsl_i = 0; curr_t <= t; ++tsl_i)
		curr_t += timeslice_sizes_pattern[tsl_i % tsl_num];
	return curr_t;
}

bool is_timeslice_boundary(int t, const timeslice_sizes& timeslice_sizes_pattern) {
	int curr_t = 0;
	while (curr_t <= t) {
		if (t == curr_t)
			return true;

		curr_t = next_timeslice_boundary(curr_t, timeslice_sizes_pattern);
	}
	return false;
}

/**
 * @param operator_t_extent					size of the operator in temporal direction
 * @param boundary_forbidden_time_coords	for each set in no_boundary_times, at least of the times in that set must
 * 											not be a timeslice boundary (i.e. elements in the outer set are related by
 * 											&&, elements in the inner sets by ||)
 * @return	set of time coordinates where an operator can be placed to start and end on a timeslice
 * 			boundary, and fulfill the required conditions given by the parameters
 */
auto find_possible_operator_time_coords(t_extent operator_t_extent,
		const timeslice_sizes& timeslice_sizes_pattern, const std::set<std::set<t_coord> >& boundary_forbidden_time_coords) {

	std::set<t_coord> possible_times;

	for (int operator_pos_t = 0; operator_pos_t < std::accumulate(timeslice_sizes_pattern.begin(), timeslice_sizes_pattern.end(), 0);
			operator_pos_t = next_timeslice_boundary(operator_pos_t, timeslice_sizes_pattern)) {

		bool non_boundaries_fulfilled = true;
		for (const auto& one_of_t_not_a_boundary : boundary_forbidden_time_coords) {
			bool found_non_boundary = false;
			for (const int t : one_of_t_not_a_boundary) {
				if (!is_timeslice_boundary(operator_pos_t + t, timeslice_sizes_pattern)) {
					found_non_boundary = true;
					break;
				}
			}
			if (!found_non_boundary) {
				non_boundaries_fulfilled = false;
				break;
			}
		}

		if (non_boundaries_fulfilled && is_timeslice_boundary(operator_pos_t + operator_t_extent, timeslice_sizes_pattern))
			possible_times.insert(operator_pos_t);
	}

	return possible_times;
}

auto try_place_operator(const timeslice_sizes& timeslice_sizes_pattern,
		const std::set<t_extent>& t_extents, observable_mode mode) {

	timeslice_setups::mapped_type possible_extents;

	for (const t_extent operator_t_extent : t_extents) {
		std::set<t_coord> possible_operator_t_coords = find_possible_operator_time_coords(operator_t_extent,
				timeslice_sizes_pattern, get_boundary_forbidden_time_coords(operator_t_extent, mode));

		if (!possible_operator_t_coords.empty())
			possible_extents[operator_t_extent] = possible_operator_t_coords;
	}

	return possible_extents;
}

void scan_operator_placements(const timeslice_sizes& timeslice_sizes_pattern,
		const std::set<t_extent>& t_extents, observable_mode mode, const std::set<int>& total_pattern_sizes,
		timeslice_setups& setups_result) {

	if (total_pattern_sizes.empty() || total_pattern_sizes.count(
			std::accumulate(timeslice_sizes_pattern.begin(), timeslice_sizes_pattern.end(), 0))) {
		auto possible_extents = try_place_operator(timeslice_sizes_pattern, t_extents, mode);
		if (!possible_extents.empty())
			setups_result[timeslice_sizes_pattern] = possible_extents;
	}
}

timeslice_sizes cycle_sizes(const timeslice_sizes& sizes) {

	timeslice_sizes cycled_sizes(sizes);

	if (sizes.size() < 2)
		return sizes;

	cycled_sizes.insert(cycled_sizes.begin(), cycled_sizes.back());
	cycled_sizes.pop_back();

	return cycled_sizes;
}

void remove_equal_patterns(timeslice_setups& setups) {
	for (const auto& [sizes, extents] : setups) {
		(void) extents;
		timeslice_sizes cycled = sizes;
		do {
			timeslice_sizes reverse = cycled;
			std::reverse(reverse.begin(), reverse.end());
			if (reverse != sizes)
				if (auto rev_pos = setups.find(reverse); rev_pos != setups.end())
					setups.erase(rev_pos);

			cycled = cycle_sizes(cycled);
			if (cycled != sizes)
				if (auto cyc_pos = setups.find(cycled); cyc_pos != setups.end())
					setups.erase(cyc_pos);

		} while (cycled != sizes);
	}
}

}

timeslice_setups find_sizes(const std::vector<std::pair<int, int> >& size_limits, const std::set<int>& total_pattern_sizes,
		const std::set<t_extent>& t_extents, observable_mode mode) {

	timeslice_setups setups;

	tools::helper::nest_for(size_limits, scan_operator_placements,
			t_extents, mode, total_pattern_sizes, setups);

	remove_equal_patterns(setups);

	return setups;
}

void filter(timeslice_setups& setups, const std::set<t_extent>& required_extents) {

	for (auto it = setups.begin(); it != setups.end();) {
		bool has_required = true;
		for (t_extent ext : required_extents)
			if (!it->second.count(ext)) {
				has_required = false;
				break;
			}
		if (has_required)
			++it;
		else
			it = setups.erase(it);
	}
}

void print(const timeslice_setups& setups) {
	using tools::io_tools::operator<<;

	for (const auto& pattern : setups) {
		std::cout << "[" << std::accumulate(pattern.first.begin(), pattern.first.end(), 0) << "] " << pattern.first << " : ";
		for (const auto& t_extent : pattern.second)
			std::cout << t_extent.first << " " << t_extent.second << " ";
		std::cout << "\n";
	}
}

}
