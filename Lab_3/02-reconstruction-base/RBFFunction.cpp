#include <iostream>
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "RBFFunction.h"
#include <typeinfo>

using namespace Eigen;

/* Initialize everything to be able to compute the implicit distance to the reconstructed
   point cloud at arbitrary points that are close enough to the point cloud. As should be
   obvious by the name of the class, the distance has to be computed using RBFs.
 */

void RBFFunction::init(const PointCloud *pointCloud, float standardDeviation, float supportRadius)
{
	cloud = pointCloud;
	c_std = standardDeviation;
	supp_r = supportRadius;
	is_success = calc_c();
}

bool RBFFunction::calc_c(){
	size_t c_size = cloud->getPoints().size();
	SparseMatrix<double> A(3*c_size, 3*c_size);
	VectorXd v(3*c_size);
	float d = 0.1;
	
	vector<glm::vec3> cloud_full(3*c_size);
	vector<glm::vec3> n = cloud->getNormals();
	for (int i = 0; i < c_size; i++)
	{
		glm::vec3 Pi = cloud->getPoints()[i];
		cloud_full[3*i] = Pi;            v[3*i] = 0;
		cloud_full[3*i+1] = Pi + d*n[i]; v[3*i+1] = d;
		cloud_full[3*i+2] = Pi - d*n[i]; v[3*i+2] = -d;
	}
	
	for (int i = 0; i < 3*c_size; i++)
	{
		glm::vec3 Pi = cloud_full[i];
		for (int j = 0; j <= i; j++)
		{
			glm::vec3 Pj = cloud_full[j];
			float ij_dist = glm::distance(Pi, Pj);
			if (ij_dist < supp_r)//3*c_std)
			{
				double RBF = exp(pow(ij_dist,2) / (2*c_std));
				A.coeffRef(i,j) = RBF;
				A.coeffRef(j,i) = RBF;
			}
		}
		std::cout << i << std::endl;
	}
	
	BiCGSTAB<SparseMatrix<double> > solver;
	solver.compute(A);
	auto c_temp = solver.solve(v);
	std::cout << (solver.info()!=Success) << std::endl;
	if(solver.info()!=Success)
		return false;
	
	const VectorXd c_temp_2 = c_temp;
	
	c.resize(c_size);
	std::cout << "i" << std::endl;
	std::cout << typeid(c_temp_2).name() << std::endl;
	for (int i = 0; i < c_size; i++)
	{
		//std::cout << c_temp[3*i] << std::endl;
		c[i] = c_temp_2(3*i);
		//std::cout << i << " : " << c[i] << std::endl;
		//std::cout << i << std::endl;
	}
	std::cout << c[0] << std::endl;
	
	return true;
}

/* This operator returns a boolean that if true signals that the value parameter
   has been modified to contain the value of the RBF implicit distance at point P.
 */

bool RBFFunction::operator()(const glm::vec3 &P, float &value) const
{
	if (!is_success)
		return false;
	
	size_t c_size = cloud->getPoints().size();
	
	value = 0;
	for (int i = 0; i < c_size; i++)
	{
		glm::vec3 Pi = cloud->getPoints()[i];
		float ij_dist = glm::distance(P, Pi);
		//std::cout << "yy : " << c[i] << std::endl;
		value += exp(pow(ij_dist,2) / (2*c_std)) * c[i];
		//std::cout << i << std::endl;
	}
	
	//std::cout << "d : " << value << std::endl;
	
	return true;
}




