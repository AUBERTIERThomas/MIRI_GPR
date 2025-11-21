#include "SimpleDistance.h"


/* Initialize everything to be able to compute the implicit distance of [Hoppe92] 
   at arbitrary points that are close enough to the point cloud.
 */

void SimpleDistance::init(const PointCloud *pointCloud, float samplingRadius)
{
	cloud = pointCloud;
	knn.setPoints(&pointCloud->getPoints());
	radius = samplingRadius;
}


/* This operator returns a boolean that if true signals that the value parameter
   has been modified to contain the value of the implicit function of [Hoppe92]
   at point P.
 */

bool SimpleDistance::operator()(const glm::vec3 &P, float &value) const
{
	vector<size_t> neighbors;
	vector<float> dists_squared;
	knn.getKNearestNeighbors(P, 1, neighbors, dists_squared);
	
	glm::vec3 Pi = cloud->getPoints()[neighbors[0]];
	glm::vec3 ni = cloud->getNormals()[neighbors[0]];
	float temp = glm::dot(P - Pi, ni);
	glm::vec3 z = Pi - temp*ni;
	if (glm::distance(z, Pi) <= radius)
	{
		value = temp;
		return true;
	}

	return false;
}






