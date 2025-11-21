#include <iostream>
#include <algorithm>
#include "IterativeClosestPoint.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace Eigen;

void IterativeClosestPoint::setClouds(PointCloud *pointCloud1, PointCloud *pointCloud2)
{
	cloud1 = pointCloud1;
	cloud2 = pointCloud2;
	knn.setPoints(&(cloud1->getPoints()));
	isBorder.assign(cloud1->getPoints().size(), false);
	colors1 = cloud1->getColors();
	correspondence = new vector<int>(cloud2->getPoints().size());
}

// This method should mark the border points in cloud 1. It also changes their color (for example to red).
// You will need to add an attribute to this class that stores this property for all points in cloud 1. 

void IterativeClosestPoint::markBorderPoints()
{
	size_t c1_l = cloud1->getPoints().size();
	unsigned int K = 20;
	float angleThreshold = 3.141592654;
	int cpt = 0;
	for (int point = 0; point < c1_l; point++)
	{
		vector<size_t> neighbors;
		vector<float> dists_squared;
		knn.getKNearestNeighbors(cloud1->getPoints()[point], K, neighbors, dists_squared);
		
		glm::vec3 centroid(0.0, 0.0, 0.0);
		size_t l = neighbors.size();
		for (int i = 0; i < l; i++)
			centroid += cloud1->getPoints()[neighbors[i]];
		centroid /= l;
		
		vector<glm::vec3> centered_points(l);
		for (int i = 0; i < l; i++)
			centered_points[i] = cloud1->getPoints()[neighbors[i]] - centroid;
		
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
		
		auto eigv = eigensolver.eigenvectors();
		auto col1 = eigv.col(0);
		glm::vec3 v1; v1.x = col1(0); v1.y = col1(1); v1.z = col1(2);
		auto col2 = eigv.col(1);
		glm::vec3 v2; v2.x = col2(0); v2.y = col2(1); v2.z = col2(2);
		auto col3 = eigv.col(2);
		glm::vec3 v3; v3.x = col3(0); v3.y = col3(1); v3.z = col3(2);
		
		vector<glm::vec3> nn_frame(K);
		float angles[K];
		for (int i = 0; i < K; i++)
		{
			auto pi_tilde = cloud1->getPoints()[neighbors[i]]-cloud1->getPoints()[point];
			glm::vec3 nn_t (glm::dot(pi_tilde, v1), glm::dot(pi_tilde, v2),glm::dot(pi_tilde, v3));
			
			nn_frame[i] = nn_t;
			angles[i] = atan2(nn_t[1], nn_t[0]);
		}
		std::sort(angles, angles+K);
		
		float maxDeltaAlpha = 0.0;
		for (int i = 0; i < K-1; i++)
		{
			float deltaAlpha = angles[i+1] - angles[i];
			if (maxDeltaAlpha < deltaAlpha) maxDeltaAlpha = deltaAlpha;
		}
		
		float deltaAlpha = angles[0] - angles[K-1] + 6.283185307;
		/*for (int i = 0; i < K; i++)
			std::cout << angles[i] << ",";
		std::cout << angles[0] << " | " << angles[K-1] << " | " << deltaAlpha << std::endl;*/
		if (maxDeltaAlpha < deltaAlpha) maxDeltaAlpha = deltaAlpha;
		
		if (maxDeltaAlpha >= angleThreshold)
		{
			isBorder[point] = true; cpt++;
			colors1[point].r = 1.0; colors1[point].g = 0.0; colors1[point].b = 0.0; colors1[point].a = 1.0; //RED
		}
	}
	std::cout << ">>> " << cpt << std::endl;
}


// This method should compute the closest point in cloud 1 for all non border points in cloud 2. 
// This correspondence will be useful to compute the ICP step matrix that will get cloud 2 closer to cloud 1.
// Store the correspondence in this class as the following method is going to need it.
// As it is evident in its signature this method also returns the correspondence. The application draws this if available.

vector<int> *IterativeClosestPoint::computeCorrespondence()
{
	size_t c1_s = cloud1->getPoints().size();
	size_t c2_s = cloud2->getPoints().size();
	//vector<int>* corr = new vector<int>();//(c2_s);
	vector<int>& corr_values = *correspondence;
	for (int point = 0; point < c2_s; point++)
	{
		int minInd = 0;
		float minNorm = glm::distance(cloud2->getPoints()[point],cloud1->getPoints()[0]);
		for (int i = 1; i < c1_s; i++)
		{
			float norm = glm::distance(cloud2->getPoints()[point],cloud1->getPoints()[i]);
			if ((minNorm > norm) && (!isBorder[i]))
			{
				minNorm = norm;
				minInd = i;
			}
		}
		//corr->push_back(minInd);
		corr_values[point] = minInd;
		//std::cout << point << " | " << minInd << " | " << minNorm << std::endl;
	}
	//correspondence = corr;
	return correspondence;//return NULL;
}


// This method should compute the rotation and translation of an ICP step from the correspondence
// information between clouds 1 and 2. Both should be encoded in the returned 4x4 matrix.
// To do this use the SVD algorithm in Eigen.

glm::mat4 IterativeClosestPoint::computeICPStep()
{
	std::cout << "a" << std::endl;
	vector<int>& corr = *correspondence;
	std::cout << "a'" << std::endl;
	glm::vec3 centroid1(0.0, 0.0, 0.0);
	glm::vec3 centroid2(0.0, 0.0, 0.0);
	size_t c2_s = cloud2->getPoints().size();
	vector<glm::vec3> centered_points_1(c2_s);
	vector<glm::vec3> centered_points_2(c2_s);
	//size_t c1_s = cloud1->getPoints().size();
	std::cout << "a''" << std::endl;
	for (int i = 0; i < c2_s; i++){
		//std::cout << i << " | " << corr[i] << std::endl;
		centroid1 += cloud1->getPoints()[corr[i]];}
	centroid1 /= c2_s;
	std::cout << "b" << std::endl;
	
	for (int i = 0; i < c2_s; i++)
		centered_points_1[corr[i]] = cloud1->getPoints()[corr[i]] - centroid1;
	std::cout << "c" << std::endl;
	
	for (int i = 0; i < c2_s; i++)
		centroid2 += cloud2->getPoints()[i];
	centroid2 /= c2_s;
	std::cout << "d" << std::endl;
	
	for (int i = 0; i < c2_s; i++){
		//std::cout << i << std::endl;
		centered_points_2[i] = cloud2->getPoints()[i] - centroid2;}
	std::cout << "e" << std::endl;
	MatrixXd Pt(c2_s, 3);
	for (int r = 0; r < c2_s; r++)
		for (int c = 0; c < 3; c++)
			Pt(r,c) = centered_points_1[r][c];
	std::cout << "f" << std::endl;
	MatrixXd Q(3, c2_s);
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < c2_s; c++)
			Q(r,c) = centered_points_2[c][r];
	std::cout << "g" << std::endl;
	MatrixXd S = Q*Pt;
	std::cout << "h" << std::endl;
	JacobiSVD<MatrixXd> svd; svd.compute(S, ComputeThinU | ComputeThinV);
	auto U = svd.matrixU();
	auto V = svd.matrixV();
	std::cout << "i" << std::endl;
	auto R = V*U.transpose();
	std::cout << "j" << std::endl;
	Vector3d cen1(centroid1.x, centroid1.y, centroid1.z);
	Vector3d cen2(centroid2.x, centroid2.y, centroid2.z);
	Vector3d trans(R(0,0)*cen2(0)+R(0,1)*cen2(1)+R(0,2)*cen2(2), R(1,0)*cen2(0)+R(1,1)*cen2(1)+R(1,2)*cen2(2), R(2,0)*cen2(0)+R(2,1)*cen2(1)+R(2,2)*cen2(2));
	Vector3d t = cen1 - trans;
	std::cout << "k" << std::endl;
	glm::mat4 res(R(0,0),R(1,0),R(2,0),0.0, R(0,1),R(1,1),R(2,1),0.0, R(0,2),R(1,2),R(2,2),0.0, t(0),t(1),t(2),0.0);
	std::cout << "l" << std::endl;
	//std::cout << R(0,0) << " " << R(1,0) << " " << R(2,0) << " " << 0.0 << " " << R(0,1) << " " << R(1,1) << " " << R(2,1) << " " << 0.0 << " " << R(0,2) << " " << R(1,2) << " " << R(2,2) << " " << 0.0 << " " << t(0) << " " << t(1) << " " << t(2) << " " << 0.0 << std::endl;
	std::cout << glm::to_string(res) << std::endl;
	return res;
	//glm::mat4 haha(1.0,0,0,0,1.0,0,0,0,0,0,1.0,0,0,0,0,0); return haha;
}


// This method should perform the whole ICP algorithm with as many steps as needed.
// It should stop when maxSteps are performed, when the Frobenius norm of the transformation matrix of
// a step is smaller than a small threshold, or when the correspondence does not change from the 
// previous step.

vector<int> *IterativeClosestPoint::computeFullICP(unsigned int maxSteps)
{
	std::cout << "hello" << std::endl;
	float epsilon = 0.01;
	size_t c2_s = cloud2->getPoints().size();
	int cpt = 0;
	double F_norm = 1; //Arbitrary initial value
	std::cout << "bjr" << std::endl;
	
	do{
		correspondence = computeCorrespondence();
	
		glm::mat4 res = computeICPStep();
		glm::mat4 I(1.0);
		glm::mat4 res_I = res - I;
		F_norm = 0;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				F_norm += res_I[i][j]*res_I[i][j];
		
		for (int i = 0; i < c2_s; i++)
		{
			glm::vec3 v2 = cloud2->getPoints()[i];
			glm::vec3 new_v2(v2[0]*res[0][0]+v2[1]*res[0][1]+v2[2]*res[0][2], 
							 v2[0]*res[1][0]+v2[1]*res[1][1]+v2[2]*res[1][2], 
							 v2[0]*res[2][0]+v2[1]*res[2][1]+v2[2]*res[2][2]);
			cloud2->getPoints()[i] = new_v2;
		}
		cpt++;
		std::cout << cpt << " ";
		
	} while ((F_norm > epsilon) && (cpt < maxSteps));
	
	return correspondence;//return NULL;
}





