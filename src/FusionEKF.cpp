#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // H matrix - laser
  H_laser_ << 1, 0, 0, 0,
      0, 1, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  tools = Tools();
  initKalman();

}

/*
 * initializer
 */

void FusionEKF::initKalman() {

  //create a 4D state vector, we don't know yet the values of the x state
  VectorXd x_ = VectorXd(4);

  //state covariance matrix P
  MatrixXd P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 10000, 0,
      0, 0, 0, 10000;

  //measurement covariance
  MatrixXd R_ = MatrixXd(2, 2);
  R_ << 0.0225, 0,
      0, 0.0225;

  //measurement matrix
  MatrixXd H_ = MatrixXd(2, 4);
  H_ << 1, 0, 0, 0,
      0, 1, 0, 0;

  //the initial transition matrix F_
  MatrixXd F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
      0, 1, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1;

  // initial process covariance Matrix Q
  MatrixXd Q_ = MatrixXd(4, 4);
  Q_ << 0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0;

  ekf_ = KalmanFilter();
  ekf_.Init(x_, P_, F_, H_, R_, Q_);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    // record time
    previous_timestamp_ = measurement_pack.timestamp_;

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    float m0 = measurement_pack.raw_measurements_[0];
    float m1 = measurement_pack.raw_measurements_[1];

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      ekf_.x_ << m0*cos(m1), m0*sin(m1), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << m0, m1, 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //1. Modify the F matrix so that the time is integrated
  ekf_.F_ << 1, 0, dt, 0,
      0, 1, 0, dt,
      0, 0, 1, 0,
      0, 0, 0, 1;

  //2. Set the process covariance matrix Q
  float dt4 = 0.25*pow(dt, 4);
  float dt3 = 0.5*pow(dt, 3);
  float dt2 = pow(dt, 2);
  float ax = 9;
  float ay = 9;

  ekf_.Q_ <<
      dt4*ax, 0, dt3*ax, 0,
      0, dt4*ay, 0, dt3*ay,
      dt3*ax, 0, dt2*ax, 0,
      0, dt3*ay, 0, dt2*ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.Update(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

}
