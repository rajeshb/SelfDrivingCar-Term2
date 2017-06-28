#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
    
    ///* initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;
    
    ///* if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;
    
    ///* if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;
    
    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;
    
    ///* state covariance matrix
    MatrixXd P_;
    
    ///* augmented sigma points matrix
    MatrixXd Xsig_aug_;
    
    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;
    
    ///* State dimension
    int n_x_;
    
    ///* Augmented state dimension
    int n_aug_;
    
    ///* Sigma point spreading parameter
    double lambda_;
    
    ///* S matrix
    MatrixXd S_;
    
    ///* Z sigma
    MatrixXd Zsig_;
    
    ///* prediction state
    VectorXd z_pred_;
    
    ///* NIS for Radar
    double nis_radar_;
    
    ///* NIS for Laser
    double nis_laser_;
    
    ///* time when the state is true, in us
    long long time_us_;
    
    ///* Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;
    
    ///* Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;
    
    ///* Laser measurement noise standard deviation position1 in m
    double std_laspx_;
    
    ///* Laser measurement noise standard deviation position2 in m
    double std_laspy_;
    
    ///* Radar measurement noise standard deviation radius in m
    double std_radr_;
    
    ///* Radar measurement noise standard deviation angle in rad
    double std_radphi_;
    
    ///* Radar measurement noise standard deviation radius change in m/s
    double std_radrd_ ;
    
    ///* Weights of sigma points
    VectorXd weights_;
    
    /**
     * Constructor
     */
    UKF();
    
    /**
     * Destructor
     */
    virtual ~UKF();
    
    /**
     * ProcessMeasurement
     * @param meas_package The latest measurement data of either radar or laser
     */
    void ProcessMeasurement(MeasurementPackage meas_package);
    
    /**
     * Prediction Predicts sigma points, the state, and the state covariance
     * matrix
     * @param delta_t Time between k and k+1 in s
     */
    void Prediction(double delta_t);
    
    /**
     * Updates the state and the state covariance matrix using a laser measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateLidar(MeasurementPackage meas_package);
    
    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateRadar(MeasurementPackage meas_package);
    
protected:
    /**
     * Generate augmented sigma points
     */
    void AugmentedSigmaPoints();
    
    /**
     * Sigma point prediction
     * @param double delta_t
     */
    void SigmaPointPrediction(double delta_t);
    
    /**
     * Mean and Covariance matrix of predicted state
     */
    void PredictMeanAndCovariance();
    
    /**
     * Predict Measurement
     * @param int n_z == 3 for radar, n_z == 2 for laser
     */
    void PredictMeasurement(int n_z);
    
    /**
     * Update State
     * @param VectorXd z
     * @param int n_z == 3 for radar, n_z == 2 for laser
     */
    void UpdateState(VectorXd z, int n_z);
};

#endif /* UKF_H */
