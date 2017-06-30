#include "PID.h"

/*
* PID class.
*/

PID::PID() {p_error_initialized = false; }

PID::~PID() {}

void PID::Init(double Kp, double Ki, double Kd) {
    this->Kp = Kp;
    this->Ki = Ki;
    this->Kd = Kd;
    
    p_error = 0.0;
    i_error = 0.0;
    d_error = 0.0;
}

void PID::UpdateError(double cte) {
    if (!p_error_initialized) {
        this->p_error = cte;
        p_error_initialized = true;
    }
    this->d_error = cte - this->p_error;
    this->i_error += cte;
    this->p_error = cte;
}

double PID::TotalError() {
    return -Kp * p_error - Kd * d_error - Ki * i_error;
}

