/*
 * multilevel.cc
 *
 *  Created on: 19.07.2019
 *      Author: reisinger
 */

#include <map>
#include <vector>

using sublattice_operator = double;

//example : compute <[T(0)T(1)][[T(2)][T(3)]]>

std::map<int, std::vector<int> > sublattice_partitions {
		{ 0, { 1, 2 } },
		{ 2, { 21, 22 } }
};

std::map<int, int> sublattice_measurement_num {
		{ 0, 20 },
		{ 1, 20 },
		{ 2, 20 },
		{ 21, 20 },
		{ 22, 20 }
};

std::map<int, std::vector<int> > sublattice_times {
		{ 1, { 0, 1 } },
		{ 21, { 2 } },
		{ 22, { 3 } }
};

sublattice_operator* compute_sublattice_operator(int sublattice_key, sublattice_operator* workspace) {
	if(sublattice_partitions.count(sublattice_key)) {

	}
}
