#include "NormalEstimator.h"
#include "NearestNeighbors.h"
#include "eigen/Eigen/Dense"
#include <iostream>

using namespace Eigen;

// This method has to compute a normal per point in the 'points' vector and put it in the 
// 'normals' vector. The 'normals' vector already has the same size as the 'points' vector. 
// There is no need to push_back new elements, only to overwrite them ('normals[i] = ...;')
// In order to compute the normals you should use PCA. The provided 'NearestNeighbors' class
// wraps the nanoflann library that computes K-nearest neighbors effciently. 

void NormalEstimator::computePointCloudNormals(const vector<glm::vec3> &points, vector<glm::vec3> &normals)
{

	float l_full = points.size();
	unsigned int K = 20;
	for (int point = 0; point < l_full; point++)
	{
		vector<size_t> neighbors;
		vector<float> dists_squared;
		NearestNeighbors nn;
		nn.getKNearestNeighbors(points[point], K, neighbors, dists_squared);
		
		glm::vec3 centroid(0.0, 0.0, 0.0);
		float l = neighbors.size();
		for (int i = 0; i < l; i++)
			centroid += points[neighbors[i]];
		centroid /= l;
		
		vector<glm::vec3> centered_points(l);
		for (int i = 0; i < l; i++)
			centered_points[i] = points[neighbors[i]] - centroid;
		
		Matrix3f C;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c <= r; c++)
			{
				float v = 0.0;
				for (int p = 0; p < l; p++)
					v += centered_points[p][r]*centered_points[p][c];
				C(r,c) = v;
				if (r != c)
					C(c,r) = v;
			}
		
		SelfAdjointEigenSolver<Matrix3f> eigensolver(C);
		if (eigensolver.info() != Success)
		{
			std::cout << "Can't compute eigenvalues !" << std::endl;
			abort();
		}
		
		//vector<float> eig_values = eigensolver.eigenvalues().col(0);
		//stable_sort(eig_values.begin(), eig_values.end());
		//vector<float> eig_values_2 = eigensolver.eigenvalues().col(0);
		//auto it = find(eig_values_2.begin(), eig_values_2.end(), eig_v[0]);
		auto eigv_matrix = eigensolver.eigenvectors();
		glm::vec3 norm_v;
		norm_v.x = eigv_matrix(0,2);
		norm_v.y = eigv_matrix(1,2);
		norm_v.z = eigv_matrix(2,2);
		normals[point] = norm_v;
	}
}


