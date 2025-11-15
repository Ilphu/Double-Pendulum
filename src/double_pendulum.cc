/**
 * @file double_pendulum.cc
 * @author Garrett Rhoads
 * @brief 
 */

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <stdexcept>
#include "double_pendulum.h"

/**
 * @brief default contructor
 */
DoublePendulum::DoublePendulum() {}

/**
 * @brief creates a double pendulum object and initializes all attributes
 */
DoublePendulum::DoublePendulum(const double& theta1,  const double& theta2, 
                               const double& mass1,   const double& mass2, 
                               const double& length1, const double& length2) {
    _mass1 = mass1;
    _length1 = length1;

    _mass2 = mass2;
    _length2 = length2;

    _theta1 = theta1;
    _theta2 = theta2;
    _omega1 = 0.0;
    _omega2 = 0.0;
    _d_omega1 = 0.0;
    _d_omega2 = 0.0;

    _state.push_back(_theta1);
    _state.push_back(_theta2);
    _state.push_back(_omega1);
    _state.push_back(_omega2);
}

DoublePendulum::~DoublePendulum() {}

void DoublePendulum::init(const double& theta1,  const double& theta2, 
                          const double& mass1,   const double& mass2, 
                          const double& length1, const double& length2) {
    
    _mass1 = mass1;
    _length1 = length1;

    _mass2 = mass2;
    _length2 = length2;

    _theta1 = theta1;
    _theta2 = theta2;
    _omega1 = 0.0;
    _omega2 = 0.0;
    _d_omega1 = 0.0;
    _d_omega2 = 0.0;

    _state.push_back(_theta1);
    _state.push_back(_theta2);
    _state.push_back(_omega1);
    _state.push_back(_omega2);
}

double DoublePendulum::get_mass1() const { return _mass1; }
double DoublePendulum::get_length1() const { return _length1; }
double DoublePendulum::get_mass2() const { return _mass2; }
double DoublePendulum::get_length2() const { return _length2; }
vector<double> DoublePendulum::get_state() const { return _state; }
void DoublePendulum::set_mass1(const double& new_mass1) { _mass1 = new_mass1; }
void DoublePendulum::set_length1(const double& new_length1) { _length1 = new_length1; }
void DoublePendulum::set_mass2(const double& new_mass2) { _mass2 = new_mass2; }
void DoublePendulum::set_length2(const double& new_length2) { _length2 = new_length2; }
void DoublePendulum::set_state(const vector<double>& new_state) { _state[0] = new_state[0]; _state[1] = new_state[1];  
                                                                 _state[2] = new_state[2]; _state[3] = new_state[3]; }

vector<double> DoublePendulum::diff_eq(const vector<double>& y, double t) {
    t++; // to make cmake happy
    double theta1 = y[0];
    double theta2 = y[1];
    double omega1 = y[2];
    double omega2 = y[3];
    double num1;
    double num2;
    double num3;
    double num4;
    double den;

    double f1 = omega1;
    double f2 = omega2;

    num1 = -g * ((2 * _mass1) + _mass2) * sin(theta1);
    num2 = -_mass2 * g * sin(theta1 - (2 * theta2));
    num3 = -2 * sin(theta1 - theta2) * _mass2;
    num4 = ((omega2 * omega2 * _length2) + (omega1 * omega1 * _length1 * cos(theta1 - theta2)));
    den  = (_length1 * ((2 * _mass1) + _mass2 - (_mass2 * cos((2 * theta1) - (2 * theta2)))));
    double f3 = (num1 + num2 + (num3 * num4)) / den;

    num1 = 2 * sin(theta1 - theta2);
    num2 = (omega1 * omega1 * _length1 * (_mass1 + _mass2));
    num3 = (g * (_mass1 + _mass2) * cos(theta1));
    num4 = omega2 * omega2 * _length1 * _mass2 * cos(theta1 - theta2);
    den  = (_length2 * ((2 * _mass1) + _mass2 - (_mass2 * cos((2 * theta1) - (2 * theta2)))));
    double f4 = (num1 * (num2 + num3 + num4)) / den;

    vector<double> new_state;
    new_state.push_back(f1);
    new_state.push_back(f2);
    new_state.push_back(f3);
    new_state.push_back(f4);
    return new_state;
}

vector<double> DoublePendulum::RK4_step(const vector<double>& y, const double& t, const double& dt) {
    vector<double> temp = y;
    vector<double> new_state;

    vector<double> k1 = diff_eq(temp, t);

    for (unsigned int i = 0; i < y.size(); i++) {
        temp[i] += k1[i] * 0.5 * dt;
    }
    vector<double> k2 = diff_eq(temp, t + (0.5 * dt));

    temp = y;
    for (unsigned int i = 0; i < y.size(); i++) {
        temp[i] += k2[i] * 0.5 * dt;
    }
    vector<double> k3 = diff_eq(temp, t + (0.5 * dt));

    temp = y;
    for (unsigned int i = 0; i < y.size(); i++) {
        temp[i] += k3[i] * dt;
    }
    vector<double> k4 = diff_eq(temp, t + dt);
    
    for (unsigned int i = 0; i < k4.size(); i++) {
        new_state.push_back(dt * (k1[i] + (2 * k2[i]) + (2 * k3[i]) + k4[i]) / 6);
    }

    return new_state;
}

double DoublePendulum::get_energy() {
    double theta1 = _state[0];
    double theta2 = _state[1];
    double omega1 = _state[2];
    double omega2 = _state[3];
    double T = 0.5 * (_mass1 + _mass2) * _length1 * _length1 * omega1 * omega1 + (_mass2 / 2) * _length2 * _length2 * 
                      omega2 * omega2 + _mass2 * _length1 * _length2 * omega1 * omega2 * cos(theta1 - theta2);
    double U = -(_mass1 + _mass2) * _length1 * g * cos(theta1) - _mass2 * _length2 * g * cos(theta2);
    return T + U;
}

void DoublePendulum::update_state(double t, double dt) {
    vector<double> new_state = RK4_step(_state, t, dt);
    for (unsigned int i = 0; i < _state.size(); i++) {
        _state[i] += new_state[i];
    }
}

void DoublePendulum::update_d_omega() {
    double num1 = -g * ((2 * _mass1) + _mass2) * sin(_theta1);
    double num2 = -_mass2 * g * sin(_theta1 - (2 * _theta2));
    double num3 = -2 * sin(_theta1 - _theta2) * _mass2;
    double num4 = ((_omega2 * _omega2 * _length2) + (_omega1 * _omega1 * _length1 * cos(_theta1 - _theta2)));
    double den  = (_length1 * ((2 * _mass1) + _mass2 - (_mass2 * cos((2 * _theta1) - (2 * _theta2)))));

    _d_omega1 = (num1 + num2 + (num3 * num4)) / den;
    
    num1 = 2 * sin(_theta1 - _theta2);
    num2 = (_omega1 * _omega1 * _length1 * (_mass1 + _mass2));
    num3 = (g * (_mass1 + _mass2) * cos(_theta1));
    num4 = _omega2 * _omega2 * _length1 * _mass2 * cos(_theta1 - _theta2);
    den  = (_length2 * ((2 * _mass1) + _mass2 - (_mass2 * cos((2 * _theta1) - (2 * _theta2)))));

    _d_omega2 = (num1 * (num2 + num3 + num4)) / den;

    _omega1 += _d_omega1;
    _omega2 += _d_omega2;
    _theta1 += _omega1;
    _theta2 += _omega2;
}

vector<RectPoint> DoublePendulum::get_position(bool state) const {
    struct RectPoint pt1;
    struct RectPoint pt2;

    if (!state) {
        pt1.x = _length1 * sin(_theta1);
        pt1.y = -_length1 * cos(_theta1);

        pt2.x = pt1.x + (_length2 * sin(_theta2));
        pt2.y = pt1.y - (_length2 * cos(_theta2));
    } else {
        pt1.x = _length1 * sin(_state[0]);
        pt1.y = -_length1 * cos(_state[0]);

        pt2.x = pt1.x + (_length2 * sin(_state[1]));
        pt2.y = pt1.y - (_length2 * cos(_state[1]));
    }
    vector<RectPoint> position;
    position.push_back(pt1);
    position.push_back(pt2);
    return position;
}
