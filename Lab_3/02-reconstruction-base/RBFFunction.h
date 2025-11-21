#ifndef _RBF_FUNCTION_INCLUDE
#define _RBF_FUNCTION_INCLUDE


#include <glm/glm.hpp>
#include "ImplicitFunction.h"
#include "PointCloud.h"
#include "NearestNeighbors.h"


class RBFFunction : public ImplicitFunction
{

public:
	void init(const PointCloud *pointCloud, float standardDeviation, float supportRadius);
	bool calc_c();

	bool operator()(const glm::vec3 &P, float &value) const;

private:
	const PointCloud *cloud;
	float c_std;
	float supp_r;
	vector<double> c;
	bool is_success;
	
};


#endif // _RBF_FUNCTION_INCLUDE


