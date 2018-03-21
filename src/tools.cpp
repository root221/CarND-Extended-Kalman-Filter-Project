

#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rsme;
	VectorXd mean;
	VectorXd s(4);
	s << 0,0,0,0;
	for(int i=0; i<estimations.size(); i++){
		VectorXd err = ground_truth[i] - estimations[i];
		err = err.array() * err.array();
		s += err;
	}
	mean = s / estimations.size();
	rsme = mean.array().sqrt();
	return rsme;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  
	// Calculate the Jacobian
	MatrixXd H(3,4);

	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float c1 = px*px + py*py;
	float c2 = sqrt(c1);
	float c3 = c1 * c2;

	H << px/c2,py/c2,0,0,
		 -py/c1,px/c1,0,0,
		 py*(vx*py-vy*px)/c3,px*(vy*px-vx*py)/c3,px/c2,py/c2;


	return H;
}
