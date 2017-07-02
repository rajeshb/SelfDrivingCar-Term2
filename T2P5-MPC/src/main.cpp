#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.rfind("}]");
    if (found_null != string::npos) {
        return "";
    } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);
    
    for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
    }
    
    for (int j = 0; j < xvals.size(); j++) {
        for (int i = 0; i < order; i++) {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }
    
    auto Q = A.householderQr();
    auto result = Q.solve(yvals);
    return result;
}

// Center of gravity needed related to psi and epsi
const double Lf = 2.67;

// Latency for predicting time at actuation
const double dt_pred = 0.1;

int main() {
    uWS::Hub h;
    
    // MPC is initialized here!
    MPC mpc;
    
    h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                       uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        string sdata = string(data).substr(0, length);
        cout << sdata << endl;
        if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
            string s = hasData(sdata);
            if (s != "") {
                auto j = json::parse(s);
                string event = j[0].get<string>();
                if (event == "telemetry") {
                    // j[1] is the data JSON object
                    vector<double> ptsx = j[1]["ptsx"];
                    vector<double> ptsy = j[1]["ptsy"];
                    double px = j[1]["x"];
                    double py = j[1]["y"];
                    double psi = j[1]["psi"];
                    double v = j[1]["speed"];
                    double delta = j[1]["steering_angle"];
                    double a = j[1]["throttle"];

                    vector<double> wp_x;
                    vector<double> wp_y;
                    
                    // transform waypoints to be from car's perspective
                    double x_rotation = cos(-psi);
                    double y_rotation = sin(-psi);
                    for (int i = 0; i < ptsx.size(); i++) {
                        double dx = ptsx[i] - px;
                        double dy = ptsy[i] - py;
                        wp_x.push_back(dx * x_rotation - dy * y_rotation);
                        wp_y.push_back(dx * y_rotation + dy * x_rotation);
                    }

                    const int STATE_VECTOR_SIZE = 6;
                    const int NUM_DEGREES = 3;
                    
                    Eigen::Map<Eigen::VectorXd> ptsx_v(&wp_x[0], STATE_VECTOR_SIZE);
                    Eigen::Map<Eigen::VectorXd> ptsy_v(&wp_y[0], STATE_VECTOR_SIZE);
                    
                    auto coeffs = polyfit(ptsx_v, ptsy_v, NUM_DEGREES);
                    double cte = polyeval(coeffs, 0);
                    double epsi = -atan(coeffs[1]);

                    // Predict state 100ms in advance (dt_pred = 0.1) using kinematic model
                    // x, y and psi are zero after transformation, for 100ms prediction
                    double px_pred = 0.0 + v * cos(0.0) * dt_pred;
                    double py_pred = 0.0 + v * sin(0.0) * dt_pred;
                    double psi_pred = 0.0 + v * -delta / Lf * dt_pred;
                    double v_pred = v + a * dt_pred;
                    double cte_pred = cte + v * sin(epsi) * dt_pred;
                    double epsi_pred = epsi + v * -delta / Lf * dt_pred;
                    
                    /*
                     * Calculate steering angle and throttle using MPC.
                     *
                     * Both are in between [-1, 1].
                     *
                     */
                    Eigen::VectorXd state(6);
                    state << px_pred, py_pred, psi_pred, v_pred, cte_pred, epsi_pred;
                    auto vars = mpc.Solve(state, coeffs);
                    
                    json msgJson;
                    msgJson["steering_angle"] = vars[0];
                    msgJson["throttle"] = vars[1];
                    
                    //Display the MPC predicted trajectory
                    vector<double> mpc_x_vals;
                    vector<double> mpc_y_vals;
                    
                    //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                    // the points in the simulator are connected by a Green line
                    for (int i = 2; i < vars.size(); i ++) {
                        if (i%2 == 0) {
                            mpc_x_vals.push_back(vars[i]);
                        }
                        else {
                            mpc_y_vals.push_back(vars[i]);
                        }
                    }
                    
                    msgJson["mpc_x"] = mpc_x_vals;
                    msgJson["mpc_y"] = mpc_y_vals;
                    
                    //Display the waypoints/reference line
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;

                    //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                    // the points in the simulator are connected by a Yellow line
                    const int NUM_POINTS = 20;
                    const int DIST_BETWEEN_POINTS = 5;
                    for (double i = 0; i < NUM_POINTS; i++){
                        next_x_vals.push_back(i * DIST_BETWEEN_POINTS);
                        next_y_vals.push_back(polyeval(coeffs, i * DIST_BETWEEN_POINTS));
                    }
                    
                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;
                    
                    auto msg = "42[\"steer\"," + msgJson.dump() + "]";
                    std::cout << msg << std::endl;
                    // Latency
                    // The purpose is to mimic real driving conditions where
                    // the car does actuate the commands instantly.
                    //
                    // Feel free to play around with this value but should be to drive
                    // around the track with 100ms latency.
                    //
                    // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
                    // SUBMITTING.
                    this_thread::sleep_for(chrono::milliseconds(100));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                }
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });
    
    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                       size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });
    
    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });
    
    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });
    
    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}