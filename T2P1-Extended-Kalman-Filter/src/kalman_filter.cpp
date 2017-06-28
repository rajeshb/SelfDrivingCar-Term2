#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::abs;

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
    /*
     * KF Prediction step
     */
    //VectorXd u = VectorXd(2);
    //u << 0, 0;
    //this->x_ = this->F_ * this->x_ + u;
    
    this->x_ = this->F_ * this->x_;
    MatrixXd Ft = this->F_.transpose();
    this->P_ = this->F_ * this->P_ * Ft + this->Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /*
     * KF Measurement update step
     */
    VectorXd y = z - this->H_ * this->x_;
    
    // Common update
    UpdateCommon(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /*
     * KF Measurement update step, jacobian of F and H instead
     */
    VectorXd hx(3);
    float px = this->x_[0];
    float py = this->x_[1];
    float vx = this->x_[2];
    float vy = this->x_[3];
    float px2 = px*px;
    float py2 = py*py;
    float px2py2_sqrt = sqrt(px2+py2);
    float atan_value = 0.0;
    float eps = 0.0001;
    
    // Be aware that atan2(0.0,0.0) is undefined
    if (fabs(px) > eps || fabs(py) > eps) {
        atan_value = atan2(py, px);
    }
    
    // avoid divide by zero with small eps
    hx << px2py2_sqrt, atan_value, (px*vx + py*vy)/fmax(eps,px2py2_sqrt);
    
    VectorXd y = z - hx;
    
    // Normalize angle - add 2π or subtract 2π until the angle is between -pi and pi
    float phi = y[1];
    while (abs(phi) > 3.14) {
        if (phi > 0) {
            phi -= 2 * 3.14;
        }
        else {
            phi += 2 * 3.14;
        }
    }
    y[1] = phi;
    
    // Common update
    UpdateCommon(y);
}

void KalmanFilter::UpdateCommon(const VectorXd &y) {
    MatrixXd Ht = this->H_.transpose();
    MatrixXd PHt = this->P_ * Ht;
    MatrixXd S = this->H_ * PHt + this->R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  PHt * Si;
    
    //new state
    this->x_ = this->x_ + (K * y);
    MatrixXd I = MatrixXd::Identity(this->x_.size(), this->x_.size());
    this->P_ = (I - K * this->H_) * this->P_;
}
