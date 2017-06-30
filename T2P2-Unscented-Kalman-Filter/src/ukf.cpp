#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;
    
    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // state dimension
    n_x_ = 5;
    
    //set augmented dimension
    n_aug_ = 7;
    
    //define spreading parameter
    lambda_ = 3 - n_aug_;

    // initial state vector
    x_ = VectorXd(n_x_);
    
    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.5;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.5;
    
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
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    if (!is_initialized_) {
        
        /* Create the covariance matrix. */
        P_ <<   1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 15, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1;
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double r = meas_package.raw_measurements_[0];
            double phi = meas_package.raw_measurements_[1];
            //double r_dot = meas_package.raw_measurements_[2]; // Not used
            x_ << r * cos(phi), r * sin(phi), 0, 0, 0;
        }
        else { //if (meas_package.sensor_type_ == MeasurementPackage::LASER)
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        }
        
        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }
    
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;   //dt - expressed in seconds
    time_us_ = meas_package.timestamp_;
    
    // if dt is 0, no need to predict/update
    if (fabs(dt) > 0.0001) {
        //predict
        Prediction(dt);
    
        //update
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
            UpdateRadar(meas_package);
        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
            UpdateLidar(meas_package);
        }
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    //augment sigma points
    AugmentedSigmaPoints();
    //sigma point prediction
    SigmaPointPrediction(delta_t);
    //mean and covariance of predicted state
    PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    
    //predict measurement for laser
    PredictMeasurement(2);
    
    //update state for laser
    VectorXd z(2);
    z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
    UpdateState(z, 2);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    
    //predict measurement for radar
    PredictMeasurement(3);
    
    //update state for radar
    VectorXd z(3);
    z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
    UpdateState(z, 3);
}

/**
 * Generate augmented sigma points
 */
void UKF::AugmentedSigmaPoints() {
    
    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
    //create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    
    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;
    
    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    
    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)           = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_)    = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }
    
    //write result
    Xsig_aug_ = Xsig_aug;
}

/**
 * Sigma point prediction
 * @param double delta_t
 */
void UKF::SigmaPointPrediction(double delta_t) {
    
    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    //predict sigma points
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug_(0,i);
        double p_y = Xsig_aug_(1,i);
        double v = Xsig_aug_(2,i);
        double yaw = Xsig_aug_(3,i);
        double yawd = Xsig_aug_(4,i);
        double nu_a = Xsig_aug_(5,i);
        double nu_yawdd = Xsig_aug_(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
        }
        else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p = v_p + nu_a * delta_t;
        
        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }
    
    //save result
    Xsig_pred_ = Xsig_pred;
}

/**
 * Mean and Covariance matrix of predicted state
 */
void UKF::PredictMeanAndCovariance() {
    
    //create vector for weights
    VectorXd weights = VectorXd(2 * n_aug_ + 1);
    
    // set weights
    double weight_0 = lambda_/(lambda_ + n_aug_);
    weights(0) = weight_0;
    for (int i=1; i < 2 * n_aug_ + 1; i++) {
        double weight = 0.5/(n_aug_ + lambda_);
        weights(i) = weight;
    }
    
    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_ + weights(i) * Xsig_pred_.col(i);
    }
    
    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P_ = P_ + weights(i) * x_diff * x_diff.transpose();
    }
}

/**
 * Predict Measurement
 * @param int n_z == 3 for radar, n_z == 2 for laser
 */
void UKF::PredictMeasurement(int n_z) {
    
    //set vector for weights
    VectorXd weights = VectorXd(2 * n_aug_ + 1);
    double weight_0 = lambda_ /(lambda_ + n_aug_);
    weights(0) = weight_0;
    for (int i=1; i< 2 * n_aug_ + 1; i++) {
        double weight = 0.5 / (n_aug_ + lambda_);
        weights(i) = weight;
    }
    
    //create matrix for sigma points in measurement space
    Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model, n_z == 3 for radar, n_z == 2 for laser
        if (n_z == 3) {
            
            // avoid divide by zero
            double r = sqrt(p_x * p_x + p_y * p_y);
            
            // atan2(0.0, 0.0) is undefined
            double phi = 0.0;
            if (fabs(p_x) > 0.0001 || fabs(p_y) > 0.0001) {
                phi = atan2(p_y, p_x);
            }
            
            Zsig_(0,i) = r;                                                     //r
            Zsig_(1,i) = phi;                                                   //phi
            Zsig_(2,i) = (p_x * v1 + p_y * v2 ) / fmax(0.0001, r);              //r_dot
        }
        else {
            Zsig_(0,i) = p_x;
            Zsig_(1,i) = p_y;
        }
    }
    
    //mean predicted measurement
    z_pred_ = VectorXd(n_z);
    z_pred_.fill(0.0);
    for (int i=0; i < 2 * n_aug_ + 1; i++) {
        z_pred_ = z_pred_ + weights(i) * Zsig_.col(i);
    }
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig_.col(i) - z_pred_;
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        S = S + weights(i) * z_diff * z_diff.transpose();
    }
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    if (n_z == 3) {
        R <<    std_radr_ * std_radr_, 0, 0,
                0, std_radphi_ * std_radphi_, 0,
                0, 0, std_radrd_* std_radrd_;
    } else if (n_z == 2) {
        R <<    std_laspx_ * std_laspx_, 0,
                0, std_laspy_ * std_laspy_;
    }
    S = S + R;
    
    //write result
    S_ = S;
}

/**
 * Update State
 * @param VectorXd z
 * @param int n_z == 3 for radar, n_z == 2 for laser
 */
void UKF::UpdateState(VectorXd z, int n_z) {
    
    //set vector for weights
    VectorXd weights = VectorXd(2 * n_aug_ + 1);
    double weight_0 = lambda_/(lambda_ + n_aug_);
    weights(0) = weight_0;
    for (int i=1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_ + lambda_);
        weights(i) = weight;
    }
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig_.col(i) - z_pred_;
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S_.inverse();
    
    //residual
    VectorXd z_diff = z - z_pred_;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    //nis calculation
    VectorXd diff = z_diff.col(0);
    VectorXd dtST = diff.transpose() * S_.inverse();
    double nis = dtST.dot(diff);
    if (n_z == 3) {
        nis_radar_ = nis;
    } else if (n_z == 2) {
        nis_laser_ = nis;
    }
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S_ * K.transpose();
}