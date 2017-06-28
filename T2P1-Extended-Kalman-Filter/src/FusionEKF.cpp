#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

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
    
    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
    0, 0.0225;
    
    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;
    
    /**
     * Finish initializing the FusionEKF.
     * Set the process and measurement noises
     */
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
    0, 1, 0, 1,
    0, 0, 1, 0,
    0, 0, 0, 1;
    
    H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;
    
    ekf_.Q_ = MatrixXd(4, 4);
    
    //set the acceleration noise components
    noise_ax = 9;
    noise_ay = 9;
    
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
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;
        
        // initial covariance matrix
        ekf_.P_ = MatrixXd(4, 4);
        ekf_.P_ <<  1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            // Polar
            float rho = measurement_pack.raw_measurements_[0];
            float phi = measurement_pack.raw_measurements_[1];
            float rho_dot = measurement_pack.raw_measurements_[2];
            // Coordinates convertion from polar to cartesian
            float x = rho * cos(phi);
            float y = rho * sin(phi);
            //float vx = rho_dot * cos(phi);
            //float vy = rho_dot * sin(phi);
            /**
             * This is something that a lot of students get it wrong.
             * phi is the direction of the object relative to our car.
             * It's not the direction in which the object is heading (i.e. heading direction).
             * So while we can perfectly calculate px and py from phi, we cannot compute vx and vy from phi.
             * We will need yaw (which is introduced in UKF) to compute vx and vy.
             * So even from radar measurement, we can only compute px and py.
             * Refer to https://goo.gl/dWj7zk
             */
            ekf_.x_ << x, y, 0 , 0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }
        
        // done initializing, no need to predict or update
        previous_timestamp_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    
    /**
     * Update the state transition matrix F according to the new elapsed time.
     - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // If dt is 0 (simultaneous measurements), you need not and should not predict again.
    if (fabs(dt) < 0.0001) {
        return;
    }
    
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    
    float dt_2 = dt * dt;
    float dt_3_2 = dt_2 * dt / 2;
    float dt_4_4 = dt_2 * dt_2 / 4;
    
    //set the process covariance matrix Q
    ekf_.Q_ <<  dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
    0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
    dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
    0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;
    
    ekf_.Predict();
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        this->Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = this->Hj_;
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }
    
    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
