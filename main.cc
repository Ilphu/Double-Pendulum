/**
 * @file main.cc
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-12-25
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h> 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))


using namespace std;

struct RectPoint {
    double x = 0;
    double y = 0;
};

struct PolarPoint {
    double r = 0;
    double theta = 0;
};

struct Pixel {
    uint8_t r = 0;
    uint8_t g = 0;
    uint8_t b = 0;
};

struct RectPoint to_rectancular(const struct PolarPoint& pp) {
    struct RectPoint pt;
    pt.x = pp.r * cos(pp.theta);
    pt.y = pp.r * sin(pp.theta);
    return pt;
}

struct PolarPoint to_polar(const struct RectPoint& pt) {
    struct PolarPoint pp;
    pp.r = sqrt((pt.x * pt.x) + (pt.y * pt.y));
    pp.theta = atan(pt.y / pt.x); // [-pi/2, pi/2]
    return pp;
}

class Graph {
public:
    Graph();
    Graph(const int& img_width, const int& img_height, const int& graph_width, const int& graph_height);
    ~Graph();

    struct Pixel get_pixel(const struct RectPoint& pt);
    void draw_pixel(const struct RectPoint& pt, const struct Pixel& px);
    void draw_line(struct RectPoint pt0, struct RectPoint pt1);

    void clear_graph();

    void write_image(const string& file_name);
private:

    int _img_width = 100;
    int _img_height = 100;
    double _graph_width = 4.0;
    double _graph_height = 4.0; 
    const int CHANNELS = 3;
    int _img_size = _img_width * _img_height * CHANNELS;;
    
    uint8_t *_img;
    uint8_t *_pix = _img;
};

Graph::Graph() {
    // _img = new uint8_t[_img_size];
}

Graph::Graph(const int& img_width, const int& img_height, const int& graph_width, const int& graph_height) {
    _img_width = img_width;
    _img_height = img_height;
    _graph_width = graph_width;
    _graph_height = graph_height;
    _img_size = _img_width * _img_height * CHANNELS;
    _img = new uint8_t[_img_size];
}

Graph::~Graph() {
    delete[] _img;
}

struct Pixel Graph::get_pixel(const struct RectPoint& pt) {
    struct Pixel px;
    int pix_x = (((pt.x / 2) + 0.5) * (_img_width - 1) * CHANNELS);
    int pix_y = ((_img_height - 1) - (((pt.y / 2) + 0.5) * (_img_height - 1)));
    _pix = _img + (((_img_width * CHANNELS * pix_y) + pix_x) - 
    ((_img_width * CHANNELS * pix_y) + pix_x) % CHANNELS);
    px.r = *(_pix + 0);
    px.g = *(_pix + 1);
    px.b = *(_pix + 2);

    return px;
}

void Graph::draw_pixel(const struct RectPoint& pt, const struct Pixel& px) {
    int pix_x = (((pt.x / _graph_width) + 0.5) * (_img_width - 1) * CHANNELS);
    int pix_y = ((_img_height - 1) - (((pt.y / _graph_height) + 0.5) * (_img_height - 1)));
    //cout << pix_x << ", " << pix_y << endl;
    _pix = _img + (((_img_width * CHANNELS * pix_y) + pix_x) - 
    ((_img_width * CHANNELS * pix_y) + pix_x) % CHANNELS);
    *(_pix + 0) = px.r;
    *(_pix + 1) = px.g;
    *(_pix + 2) = px.b;
}

void Graph::draw_line(struct RectPoint pt0, struct RectPoint pt1) {
    struct RectPoint pt;
    struct RectPoint swap;
    double dx = pt1.x - pt0.x;
    double dy = pt1.y - pt0.y;
    double m = dy / dx;

    struct Pixel px;
    px.r = 0;
    px.g = 0;
    px.b = 0;
    if (abs(m) > 1) {
        if (pt0.y > pt1.y) {
            swap.x = pt0.x;
            pt0.x = pt1.x;
            pt1.x = swap.x;
            swap.y = pt0.y;
            pt0.y = pt1.y;
            pt1.y = swap.y;
        }
        for (pt.y = pt0.y; pt.y < pt1.y; pt.y += (1.0 / _img_width)) {
            pt.x = ((pt.y - pt0.y) / m) + pt0.x;
            px.r = 0;
            px.g = 0;
            px.b = 0;
            draw_pixel(pt, px);
        }
    } else {
        if (pt0.x > pt1.x) {
            swap.x = pt0.x;
            pt0.x = pt1.x;
            pt1.x = swap.x;
            swap.y = pt0.y;
            pt0.y = pt1.y;
            pt1.y = swap.y;
        }
        for (pt.x = pt0.x; pt.x < pt1.x; pt.x += (1.0 / _img_width)) {
            pt.y = m * (pt.x - pt0.x) + pt0.y;
            px.r = 0;
            px.g = 0;
            px.b = 0;
            draw_pixel(pt, px);
        }
    }
}

void Graph::clear_graph() {
    for (_pix = _img; _pix < (_img + _img_size); _pix++) {
        *_pix = 0xFF;
    }
}

void Graph::write_image(const string& file_name) {
    stbi_write_png(file_name.c_str(), _img_width, _img_height, CHANNELS, _img, _img_width * CHANNELS);
}



class DoublePendulum {
public:
    DoublePendulum();
    DoublePendulum(const double& theta1,  const double& theta2, 
                   const double& mass1,   const double& mass2, 
                   const double& length1, const double& length2);
    ~DoublePendulum();

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

    vector<double> G(const vector<double>& y, const double& t);
    vector<double> RK4_step(const vector<double>& y, const double& t, const double& dt);
    void update_state(double t, double dt);
    void update_d_omega();
    vector<RectPoint> get_position(bool state) const;
private:
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

DoublePendulum::DoublePendulum() {}

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

double DoublePendulum::get_mass1() const {return _mass1;}
double DoublePendulum::get_length1() const {return _length1;}
double DoublePendulum::get_mass2() const {return _mass2;}
double DoublePendulum::get_length2() const {return _length2;}
vector<double> DoublePendulum::get_state() const {return _state;}
void DoublePendulum::set_mass1(const double& new_mass1) {_mass1 = new_mass1;}
void DoublePendulum::set_length1(const double& new_length1) {_length1 = new_length1;}
void DoublePendulum::set_mass2(const double& new_mass2) {_mass2 = new_mass2;}
void DoublePendulum::set_length2(const double& new_length2) {_length2 = new_length2;}
void DoublePendulum::set_state(const vector<double>& new_state) {_state[0] = new_state[0]; _state[1] = new_state[1]; 
                                                                 _state[2] = new_state[2]; _state[3] = new_state[3];}

vector<double> DoublePendulum::G(const vector<double>& y, const double& t) {
    double theta1 = y[0];
    double theta2 = y[1];
    double omega1 = y[2];
    double omega2 = y[3];
    double num1;
    double num2;
    double num3;
    double num4;
    double den;

    // double m11 = (_mass1 + _mass2) * _length1;
    // double m12 = (_mass2 * _length2 * cos(theta1 - theta2));
    // double m21 = (_length1 * cos(theta1 - theta2));
    // double m22 = (_length2);

    // double det_m = (m11 * m22) - (m12 * m21);

    // double f1 = (-1 * _mass2 * _length2 * angular_v2 * angular_v2 * sin(theta1 - theta2)) -
    //             ((_mass1 + _mass2) * g * sin(theta1));
    // double f2 = (_length1 * angular_v1 * angular_v2 * sin(theta1 - theta2)) * (g * sin(theta2));

    // double accel1 = (1 / det_m) * ((m22 * f1) - (m21 * f2));
    // double accel2 = (1 / det_m) * ((-1 * m12 * f1) + (m11 * f2));

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

    vector<double> k1 = G(temp, t);

    for (int i = 0; i < y.size(); i++) {
        temp[i] += k1[i] * 0.5 * dt;
    }
    vector<double> k2 = G(temp, t + (0.5 * dt));

    temp = y;
    for (int i = 0; i < y.size(); i++) {
        temp[i] += k2[i] * 0.5 * dt;
    }
    vector<double> k3 = G(temp, t + (0.5 * dt));

    temp = y;
    for (int i = 0; i < y.size(); i++) {
        temp[i] += k3[i] * dt;
    }
    vector<double> k4 = G(temp, t + dt);
    
    for (int i = 0; i < k4.size(); i++) {
        //cout << "here" << endl;
        new_state.push_back(dt * (k1[i] + (2 * k2[i]) + (2 * k3[i]) + k4[i]) / 6);
    }

    return new_state;
}

void DoublePendulum::update_state(double t, double dt) {
    vector<double> new_state = RK4_step(_state, t, dt);
    for (int i = 0; i < _state.size(); i++) {
        _state[i] += new_state[i];
        //cout << new_state[i] << endl;
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

int main() {
    RectPoint o;
    int img_width;
    int img_height;
    int graph_width;
    int graph_height;
    double theta1;
    double theta2;

    printf("Image resolution:\nWidth:\n");
    cin >> img_width;
    printf("Height:\n");
    cin >> img_height;
    printf("Graph size:\nWidth:\n");
    cin >> graph_width;
    printf("Height:\n");
    cin >> graph_height;
    printf("Starting angle of pedulum arms:\nPrimary:\n");
    cin >> theta1;
    printf("Secondary:\n");
    cin >> theta2;

    DoublePendulum dp(theta1, theta2, 1, 1, 1, 1);
    Graph graph(img_width, img_height, graph_width, graph_height);

    vector<RectPoint> test;

    // for (int i = 0; i < 150; i++) {
    //     test = dp.get_position();
    //     //cout << "p1 = (" << test[0].x << ", " << test[0].y << ")" << endl;
    //     //cout << "p2 = (" << test[1].x << ", " << test[1].y << ")" << endl;
    //     graph.clear_graph();
    //     graph.draw_line(o, test[0]);
    //     graph.draw_line(test[0], test[1]);
    //     string file_name = "Double_Pendulum_" + to_string(i) + ".png";
    //     graph.write_image(file_name);

    //     dp.update_d_omega();
    // }

    double dt = 0.01;
    double total_time = 100;
    double t = 0;
    int frame = 0;

    while (true) {
        test = dp.get_position(1);
        t += dt;
        dp.update_state(t, dt);
        graph.clear_graph();
        graph.draw_line(o, test[0]);
        graph.draw_line(test[0], test[1]);
        string file_name = "Double_Pendulum_" + to_string(frame) + ".png";
        graph.write_image(file_name);
        frame++;
    }
    

}