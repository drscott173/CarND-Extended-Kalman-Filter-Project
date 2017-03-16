#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

VectorXd KalmanFilter::predict_x(const VectorXd &z) {
  if (z.size() == 2) {
    return H_ * x_;
  }
  else {
    VectorXd z_pred = VectorXd(3);

    float px = x_[0];
    float py = x_[1];
    float vx = x_[2];
    float vy = x_[3];
    float h1 = sqrt(px*px + py*py);
    if ((fabs(px) < 0.0001) || (fabs(h1) < 0.0001)) {
      cout << "predict_x: Unable to predict -- zero motion" << endl;
      z_pred << nan(""), nan(""), nan("");
      return z_pred;
    }
    float h2 = atan2(py, px);
    float h3 = (px*vx + py*vy)/h1;
    if (h2 > M_PI) {
      h2 = (h2 - 2*M_PI);
    }
    z_pred << h1, h2, h3;
    return z_pred;
  }
}
    
void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = predict_x(z);
  if (isnan(z_pred[0])) return;
      
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
