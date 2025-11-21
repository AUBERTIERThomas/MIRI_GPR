#include <iostream>
#include <Eigen/Dense>
#include "MongePatch.h"
#include "glm/gtx/string_cast.hpp"

using namespace Eigen;

// Given a point P, its normal, and its closest neighbors (including itself) 
// compute a quadratic Monge patch that approximates the neighborhood of P.
// The resulting patch will be used to compute the principal curvatures of the 
// surface at point P.

void MongePatch::init(const glm::vec3 &P, const glm::vec3 &normal, const vector<glm::vec3> &closest)
{
	glm::vec3 w = -normal;
	glm::vec3 x(1,0,0);
	glm::vec3 u = glm::cross(x, w);
	glm::vec3 v = glm::cross(w, u);
	
	size_t nb_cp = closest.size();
	vector<glm::vec3> pi_list;
	
	//cout << "pi_list = " << endl;
	for(int i = 0; i < nb_cp; i++)
	{
		glm::vec3 temp(glm::dot(u, closest[i] - P), glm::dot(v, closest[i] - P), glm::dot(w, closest[i] - P));
		pi_list.push_back(temp);
		//cout << glm::to_string(temp) << endl;
	}
	MatrixXf A(6,6);
	VectorXf qi_sum(6); qi_sum << 0,0,0,0,0,0;
	//cout << "!!!" << endl;
	for (int i = 0; i < nb_cp; i++)
	{
		glm::vec3 pi = pi_list[i];
		VectorXf qi(6); qi << pi[0]*pi[0],pi[0]*pi[1],pi[1]*pi[1],pi[0],pi[1],1.;
		//cout << qi << endl;
		for (int r = 0; r < 6; r++)
			for (int c = 0; c < 6; c++)
			{
				A(c,r) = qi[c]*qi[r];
			}
		qi_sum += qi;
	}
	/*MatrixXf A(3,3);
	VectorXf qi_sum(3); qi_sum << 0,0,0;
	//cout << "!!!" << endl;
	for (int i = 0; i < nb_cp; i++)
	{
		glm::vec3 pi = pi_list[i];
		VectorXf qi(3); qi << pi[0]*pi[0],pi[0]*pi[1],pi[1]*pi[1];
		//cout << qi << endl;
		for (int r = 0; r < 3; r++)
			for (int c = 0; c < 3; c++)
			{
				A(c,r) = qi[c]*qi[r];
			}
		qi_sum += qi;
	}*/
	//cout << "..." << endl;
	//cout << "A = " << A << endl;
	//cout << qi_sum << endl;
	
	ColPivHouseholderQR<MatrixXf> dec(A);
	VectorXf s = dec.solve(qi_sum);
	cout << "s = " << s << endl;
	//cout << "|||" << endl;
	
	Matrix2f Hw; Hw << 2*s(0), s(1), s(1), 2*s(2);
	//cout << Hw << endl;
	SelfAdjointEigenSolver<Matrix2f> eigensolver(Hw);
	if (eigensolver.info() != Success)
		{
			std::cout << "Can't compute eigenvalues !" << std::endl;
			abort();
		}
		
	auto eigv = eigensolver.eigenvalues();
	K_min = eigv(0);
	K_max = eigv(1);
	//cout << "eigv = " << eigv << endl;
}

// Return the values of the two principal curvatures for this patch

void MongePatch::principalCurvatures(float &kmin, float &kmax) const
{
	kmin = K_min;
	kmax = K_max;
}


