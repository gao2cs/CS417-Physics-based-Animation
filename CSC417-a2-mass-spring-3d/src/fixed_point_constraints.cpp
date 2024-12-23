#include <fixed_point_constraints.h>
#include <algorithm>
#include <unordered_set>
// 100% Sure!

void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {

	int n = q_size / 3;
	int m = indices.size();

	P.resize(3 * (n - m), 3 * n);
	std::unordered_set<unsigned int> fixed_indcies(indices.begin(), indices.end());

	std::vector<Eigen::Triplet<double>> triplets;
	int ii = 0;
	for (int i = 0; i < q_size / 3; ++i) {
		if (fixed_indcies.find(i) != fixed_indcies.end()) {
			continue;
		}

		for (int j = 0; j < 3; ++j) {
				triplets.emplace_back(3 * ii + j, 3 * i + j, 1.0);	
		}
		ii += 1;
	}

	P.setFromTriplets(triplets.begin(), triplets.end());

}