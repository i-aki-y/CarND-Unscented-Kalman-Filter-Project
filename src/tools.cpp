#include <iostream>
#include "tools.h"

using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  size_t est_size = estimations.size();
  size_t gt_size = ground_truth.size();
  if( (est_size == 0) || (gt_size == 0)) {
    cout << "size of estimation and ground truth must be larger than 0" << endl;
    return rmse;
  }

  if(est_size != gt_size){
    cout << "sizes of estimations and ground_truth should be same" << endl;
    return rmse;
  }

  //accumulate squared residuals
  ArrayXd accum = ArrayXd(4);
  accum << 0, 0, 0, 0;
  for(int i=0; i < estimations.size(); ++i){
    VectorXd delta = ground_truth[i] - estimations[i];
    accum += delta.array() * delta.array();
  }

  //calculate the mean
  VectorXd mean = accum/estimations.size();

  //calculate the squared root
  rmse = mean.array().sqrt().matrix();

  //return the result
  return rmse;



}