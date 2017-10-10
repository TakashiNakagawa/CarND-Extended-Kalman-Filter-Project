#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace  {
  double NormalizeAngle(double angle)
  {
    return fmod(angle + M_PI, 2*M_PI) - M_PI;
  }
  
  void UpdateXP(const MatrixXd& H, const MatrixXd& R, const VectorXd& y, VectorXd& x, MatrixXd& P)
  {
    MatrixXd S = H * P * H.transpose() + R;
    MatrixXd K = P * H.transpose() * S.inverse();
    
    x = x + (K * y);
    auto x_size = x.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P = (I - K * H) * P;
  }
}// unnamed namespace

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
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;   // x = Fx + u, where here u = 0;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  // Laser uses normal Update()
  VectorXd z_pred = H_ * x_;    // Note: z is from measurement pack
  VectorXd y = z - z_pred;   

  UpdateXP(H_, R_, y, x_, P_);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // Check division by zero
  if (fabs(px*px + py*py) < 0.0001) {
    std::cout << "UpdateEKF(): CalculateJacobian() - Error - Division by Zero" << std::endl;
    return;
  }

  VectorXd z_pred = VectorXd(3);  
  float rho = sqrt(px*px + py*py); 
  float phi = atan2(py,px);         // note atan2() returns degrees in range [-pi,pi]
  float rho_dot = (px*vx + py*vy)/rho;  
  z_pred << rho, phi, rho_dot;

  // y = z_radar - h(x)_radar
  VectorXd y = z - z_pred;

  y(1) = atan2(sin(y(1)), cos(y(1)));
  y(1) = NormalizeAngle(y(1));
  
  UpdateXP(H_, R_, y, x_, P_);
}
