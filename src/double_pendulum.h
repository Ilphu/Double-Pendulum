/**
 * @file double_pendulum.h
 * @author Garrett
 * @brief 
 */

#ifndef DOUBLE_PENDULUM_H
#define DOUBLE_PENDULUM_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <stdexcept>
#include "structs.h"

using namespace std;

class DoublePendulum {
public:
    DoublePendulum();
    DoublePendulum(const double& theta1,  const double& theta2, 
                   const double& mass1,   const double& mass2, 
                   const double& length1, const double& length2);
    ~DoublePendulum();

    void init(const double& theta1,  const double& theta2, 
              const double& mass1,   const double& mass2, 
              const double& length1, const double& length2);

    double get_mass1() const;
    double get_length1() const;
    double get_mass2() const;
    double get_length2() const;
    vector<double> get_state() const;

    void set_mass1(const double& new_mass1);
    void set_length1(const double& new_length1);
    void set_mass2(const double& new_mass2);
    void set_length2(const double& new_length2);
    void set_state(const vector<double>& new_state);

    vector<double> diff_eq(const vector<double>& y, double t);
    vector<double> RK4_step(const vector<double>& y, const double& t, const double& dt);
    double get_energy();
    void update_state(double t, double dt);
    vector<RectPoint> get_position(bool state) const;
private:
    void update_d_omega();
    const double g = 9.806;

    double _mass1;
    double _length1;

    double _mass2;
    double _length2;

    double _theta1;
    double _theta2;
    double _omega1;
    double _omega2;
    double _d_omega1;
    double _d_omega2;

    vector<double> _state; // theta1 theta2, omega1 omega2,
};

#endif