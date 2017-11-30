#include "ukf.h"
//#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  //------------------------------------------------------------------------//
  //   Tunable parameter                                                    //
  //------------------------------------------------------------------------//
  //1.0, 0.5
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  //------------------------------------------------------------------------//


  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Parameters above this line are scaffolding, do not modify
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;


  P_ <<   1., 0., 0., 0., 0.,
          0., 1., 0., 0., 0.,
          0., 0., 1., 0., 0.,
          0., 0., 0., 0.5, 0.,
          0., 0., 0., 0., 0.5;


  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  //set weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(1/(2*(lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "UKF: " << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      float rho;
      float psi;

      //set the state with the initial location and zero velocity
      rho = meas_package.raw_measurements_[0];
      psi = meas_package.raw_measurements_[1];

      float px;
      float py;
      px = rho * cos(psi);
      py = rho * sin(psi);

      x_ <<   px,
              py,
              0,
              0,
              0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      float px;
      float py;
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];
      float psi;
      psi = atan2(py, px);

      x_ <<   px,
              py,
              0,
              0,
              0;

    }
    time_us_ = meas_package.timestamp_;
    //cout << "initialized." << endl;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  Prediction(dt);

  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  VectorXd x_aug = VectorXd(7);
  x_aug(0) = x_(0);
  x_aug(1) = x_(1);
  x_aug(2) = x_(2);
  x_aug(3) = x_(3);
  x_aug(4) = x_(4);
  x_aug(5) = 0;
  x_aug(6) = 0;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  //cout << "P_aug" << P_aug << endl;

  // A = sqrt(P)
  MatrixXd A = P_aug.llt().matrixL();

  //cout << "A:" << A << endl;

  //set first column of sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  Xsig_aug.col(0) = x_aug;
  //set remaining sigma points
  const double a = sqrt(lambda_ + n_aug_);
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + a * A.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - a * A.col(i);
  }


  //predict sigma points
  const double delta_t2 = delta_t * delta_t;
  //cout << "Xsig_pred " << Xsig_pred_ << endl;

  for (int i = 0; i < (1 + 2 * n_aug_); i++) {

    VectorXd x = Xsig_aug.col(i);
    const double v = x(2);
    const double psi = x(3);
    const double psi_dot = x(4);
    const double nu_a = x(5);
    const double nu_p = x(6);

    if (fabs(psi_dot) > 0.001) {
      Xsig_pred_.col(i) << x(0) + (v / psi_dot) * ( sin(psi + psi_dot * delta_t) - sin(psi)) + 0.5 * delta_t2 * cos(psi) * nu_a,
                           x(1) + (v / psi_dot) * (-cos(psi + psi_dot * delta_t) + cos(psi)) + 0.5 * delta_t2 * sin(psi) * nu_a,
                           x(2) + 0 + delta_t * nu_a,
                           x(3) + psi_dot * delta_t + 0.5 * delta_t2 * nu_p,
                           x(4) + 0 + delta_t * nu_p;
    } else {
      Xsig_pred_.col(i) << x(0) + v * cos(psi_dot) * delta_t + 0.5 * delta_t2 * cos(psi) * nu_a,
                           x(1) + v * sin(psi_dot) * delta_t + 0.5 * delta_t2 * sin(psi) * nu_a,
                           x(2) + 0 + delta_t * nu_a,
                           x(3) + 0 + 0.5 * delta_t2 * nu_p,
                           x(4) + 0 + delta_t * nu_p;

    }
    //cout << "Xsig_pred_col " << i << endl;
    //cout << Xsig_pred_.col(i) << endl;
  }

  //cout << "Xsig_pred" << Xsig_pred_ << endl;

  //predict state mean
  x_.fill(0.0);

  x_ = Xsig_pred_ * weights_;

  //predict state covariance matrix
  P_.fill(0.0);

  //cout << "x_" << x_ << endl;

  for (int i = 0; i < (2 * n_aug_ + 1); i++) {
    VectorXd dx = Xsig_pred_.col(i) - x_;

    // normalize angle
    NormalizeAngle(dx(3));

    P_ += weights_(i) * dx * dx.transpose();
  }

  //cout << "predicted" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */


  /*-----------------------------------------*/
  /*   Measurement Prediction                */
  /*-----------------------------------------*/

  const int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  z_pred.fill(0.0);
  for(int i=0; i < (2*n_aug_ + 1); i++){
    const double px = Xsig_pred_(0, i);
    const double py = Xsig_pred_(1, i);

    Zsig.col(i) << px,
                   py;

    z_pred += weights_(i) * Zsig.col(i);
  }
  //calculate mean predicted measurement
  //calculate measurement covariance matrix S
  S.fill(0.0);
  S(0, 0) = std_laspx_*std_laspx_;
  S(1, 1) = std_laspy_*std_laspy_;
  for(int i=0; i < (2*n_aug_ + 1); i++){
    VectorXd dz = Zsig.col(i)  - z_pred;

    S += weights_(i) * dz * dz.transpose();
  }

  //cout << "S:" << S << endl;

  /*-----------------------------------------*/
  /*   Update                                */
  /*-----------------------------------------*/

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i < (2*n_aug_ + 1); i++){
    VectorXd dx = Xsig_pred_.col(i) - x_;
    VectorXd dz = Zsig.col(i) - z_pred;

    // normalize
    NormalizeAngle(dx(3));

    Tc += weights_(i) * dx * dz.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //cout << "K:" << K << endl;

  // get measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],
       meas_package.raw_measurements_[1];

  //update state mean and covariance matrix
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();


  NIS_laser_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  //cout << "S:" << S << endl;

  //transform sigma points into measurement space
  z_pred.fill(0.0);
  for(int i=0; i < (2*n_aug_ + 1); i++){
    const double px =  Xsig_pred_(0, i);
    const double py =  Xsig_pred_(1, i);
    const double v  =  Xsig_pred_(2, i);
    const double psi = Xsig_pred_(3, i);

    const double rho = sqrt(px*px + py*py);
    const double phi = atan2(py, px);
    const double rho_dot = (px*cos(psi)*v + py*sin(psi)*v) / (rho);

    Zsig.col(i) << rho,
                   phi,
                   rho_dot;
    z_pred += weights_(i) * Zsig.col(i);
  }
  //calculate mean predicted measurement
  //calculate measurement covariance matrix S
  S.fill(0.0);
  S(0, 0) = std_radr_*std_radr_;
  S(1, 1) = std_radphi_*std_radphi_;
  S(2, 2) = std_radrd_*std_radrd_;
  for(int i=0; i < (2*n_aug_ + 1); i++){

    VectorXd dz = Zsig.col(i) - z_pred;

    // normalize
    NormalizeAngle(dz(1));
    S += weights_(i) * dz * dz.transpose();
  }

  /*-----------------------------------------*/
  /*   Update                                */
  /*-----------------------------------------*/

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i < (2*n_aug_ + 1); i++){
    VectorXd dx = Xsig_pred_.col(i) - x_;
    VectorXd dz = Zsig.col(i) - z_pred;

    // normalize
    NormalizeAngle(dx(3));
    NormalizeAngle(dz(1));

    Tc += weights_(i) * dx * dz.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //cout << "K:" << K << endl;

  // get measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],
       meas_package.raw_measurements_[1],
       meas_package.raw_measurements_[2];

  //update state mean and covariance matrix
  VectorXd dz = z - z_pred;
  // normalize
  NormalizeAngle(dz(1));

  x_ = x_ + K * dz;
  P_ = P_ - K * S * K.transpose();

  NIS_radar_ = dz.transpose() * S.inverse() * dz;

}

void UKF::NormalizeAngle(double& phi)
{
  phi = atan2(sin(phi), cos(phi));
}