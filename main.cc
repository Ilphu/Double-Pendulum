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
#include <thread>
#include <stdexcept>
//#include <iomapip>
#include <bits/stdc++.h> 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))


using namespace std;

const int NUM_THREADS = 8;
const int CHANNELS = 3;
const double TAU = 2 * M_PI;

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
    Graph(const int& img_width_res, const int& img_height_res, const double& graph_width, const double& graph_height);
    ~Graph();

    int get_img_width() const;
    int get_img_height() const;
    double get_graph_width() const;
    double get_graph_height() const;
    int get_img_size() const;

    // void set_img_width(const int& new_img_width);
    // void set_img_height(const int& new_img_height);
    // void set_graph_width(const double& new_graph_width);
    // void set_graph_height(const double& new_graph_height);
    
    struct Pixel get_pixel(const struct RectPoint& pt);
    struct Pixel get_pixel(const int& pix_x, const int& pix_y);
    void draw_pixel(const struct RectPoint& pt, const struct Pixel& px);
    void draw_pixel(const int& pix_x, const int& pix_y, const struct Pixel& px);
    void draw_line(struct RectPoint pt0, struct RectPoint pt1);
    void draw_line(struct RectPoint pt0, struct RectPoint pt1, const struct Pixel& px);

    void blackout();
    void whiteout();
    void bg_fill(const struct Pixel& px);

    void write_image(const string& file_name);

    Graph operator-(Graph& graph) {
        if (_img_size != graph.get_img_size()) {
            cout << "Images of different sizes" << endl;
            throw;
        }
        Graph out(_img_width, _img_height, _graph_width, _graph_height);
        struct Pixel px_this;
        struct Pixel px_graph;
        struct Pixel px_out;
        for (int i = 0; i < _img_height; i++) {
            for (int j = 0; j < _img_width; j++) {
                px_this = get_pixel(j, i);
                px_graph = graph.get_pixel(j, i);

                int r = px_this.r - px_graph.r;
                int g = px_this.g - px_graph.g;
                int b = px_this.b - px_graph.b;
            
                if (r < 0x00) {
                    r = 0x00;
                } else if (g < 0x00) {
                    g = 0x00;
                } else if (b < 0x00) {
                    b = 0x00;
                }

                px_out.r = r;
                px_out.g = g;
                px_out.b = b;
                
                out.draw_pixel(j, i, px_out);
            }
        }
        return out;
    }

    Graph& operator-=(Graph& graph) {
        if (_img_size != graph.get_img_size()) {
            cout << "Images of different sizes" << endl;
            throw;
        }
        struct Pixel px_this;
        struct Pixel px_graph;
        for (int i = 0; i < _img_height; i++) {
            for (int j = 0; j < _img_width; j++) {
                px_this = this->get_pixel(j, i);
                px_graph = graph.get_pixel(j, i);
                int r = px_this.r - px_graph.r;
                int g = px_this.g - px_graph.g;
                int b = px_this.b - px_graph.b;

                if (r < 0) { r = 0; }
                if (g < 0) { g = 0; } 
                if (b < 0) { b = 0; }

                px_this.r = r;
                px_this.g = g;
                px_this.b = b;

                this->draw_pixel(j, i, px_this);                
            }
        }
        return *this;
    }

    Graph& operator+=(Graph& graph) {
        if (_img_size != graph.get_img_size()) {
            cout << "Images of different sizes" << endl;
            throw;
        }
        struct Pixel px_this;
        struct Pixel px_graph;
        for (int i = 0; i < _img_height; i++) {
            for (int j = 0; j < _img_width; j++) {
                px_this = this->get_pixel(j, i);
                px_graph = graph.get_pixel(j, i);
                int r = px_this.r + px_graph.r;
                int g = px_this.g + px_graph.g;
                int b = px_this.b + px_graph.b;

                if (r > 0xFF) { r = 0xFF; }
                if (g > 0xFF) { g = 0xFF; } 
                if (b > 0xFF) { b = 0xFF; }

                px_this.r = r;
                px_this.g = g;
                px_this.b = b;

                this->draw_pixel(j, i, px_this);                
            }
        }
        return *this;
    }
private:

    int _img_width = 100;
    int _img_height = 100;
    double _graph_width = 4.0;
    double _graph_height = 4.0; 
    //const int CHANNELS = 3;
    int _img_size = _img_width * _img_height * CHANNELS;
    
    uint8_t *_img;
    uint8_t *_pix = _img;
};

Graph::Graph() {}

Graph::Graph(const int& img_width_res, const int& img_height_res, const double& graph_width, const double& graph_height) {
    _img_width = img_width_res;
    _img_height = img_height_res;
    _graph_width = graph_width;
    _graph_height = graph_height;
    _img_size = _img_width * _img_height * CHANNELS;
    _img = new uint8_t[_img_size];
}

Graph::~Graph() {
    delete[] _img;
}

int Graph::get_img_width() const { return _img_width; };
int Graph::get_img_height() const { return _img_height; };
double Graph::get_graph_width() const { return _graph_width; };
double Graph::get_graph_height() const { return _graph_height; };
int Graph::get_img_size() const { return _img_size; };

// void Graph::set_img_width(const int& new_img_width) {};
// void Graph::set_img_height(const int& new_img_height) {};
// void Graph::set_graph_width(const double& new_graph_width) {};
// void Graph::set_graph_height(const double& new_graph_height) {};

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

struct Pixel Graph::get_pixel(const int& pix_x, const int& pix_y) {
    struct Pixel px;

    _pix = _img + (((_img_width * CHANNELS * pix_y) + (pix_x * CHANNELS)) - 
    ((_img_width * CHANNELS * pix_y) + (pix_x * CHANNELS)) % CHANNELS);
    px.r = *(_pix + 0);
    px.g = *(_pix + 1);
    px.b = *(_pix + 2);

    return px;
}

void Graph::draw_pixel(const struct RectPoint& pt, const struct Pixel& px) {
    int pix_x = (((pt.x / _graph_width) + 0.5) * (_img_width - 1) * CHANNELS);
    int pix_y = ((_img_height - 1) - (((pt.y / _graph_height) + 0.5) * (_img_height - 1)));
    _pix = _img + (((_img_width * CHANNELS * pix_y) + pix_x) - 
    ((_img_width * CHANNELS * pix_y) + pix_x) % CHANNELS);
    *(_pix + 0) = px.r;
    *(_pix + 1) = px.g;
    *(_pix + 2) = px.b;
}

void Graph::draw_pixel(const int& pix_x, const int& pix_y, const struct Pixel& px) {
    _pix = _img + (((_img_width * CHANNELS * pix_y) + (pix_x * CHANNELS)) - 
    ((_img_width * CHANNELS * pix_y) + (pix_x * CHANNELS)) % CHANNELS);
    *(_pix + 0) = px.r;
    *(_pix + 1) = px.g;
    *(_pix + 2) = px.b;
}

void Graph::draw_line(struct RectPoint pt0, struct RectPoint pt1) {
    struct Pixel px;
    px.r = 0;
    px.g = 0;
    px.b = 0;
    draw_line(pt0, pt1, px);
}

void Graph::draw_line(struct RectPoint pt0, struct RectPoint pt1, const struct Pixel& px) {
    struct RectPoint pt;
    struct RectPoint swap;
    double dx = pt1.x - pt0.x;
    double dy = pt1.y - pt0.y;
    double m = dy / dx;
    
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
            draw_pixel(pt, px);
        }
    }
}

void Graph::blackout() {
    for (_pix = _img; _pix < (_img + _img_size); _pix++) {
        *_pix = 0x00;
    }
}

void Graph::whiteout() {
    for (_pix = _img; _pix < (_img + _img_size); _pix++) {
        *_pix = 0xFF;
    }
}

void Graph::bg_fill(const struct Pixel& px) {
    for (_pix = _img; _pix < (_img + _img_size); _pix += CHANNELS) {
        *(_pix + 0) = px.r;
        *(_pix + 1) = px.g;
        *(_pix + 2) = px.b;
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

    vector<double> G(const vector<double>& y, const double& t);
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
    for (int i = 0; i < _state.size(); i++) {
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

void render(const double& theta1, const double& theta2, const int& num_pendulums, 
            const int& res_width, const int& res_height, const int& total_frames, 
            const int& frame_rate = 60) {

    struct RectPoint origin;
    
    Graph canvas(res_width, res_height, 7.11111, 4);
    Graph trace(res_width, res_height, 7.11111, 4);
    trace.blackout();

    vector<DoublePendulum*> dp_arr;
    for (int i = 0; i < num_pendulums; i++) {
        DoublePendulum *dp = new DoublePendulum(theta1, theta2 + (0.0001 * i), 1, 1, 1, 1);
        dp_arr.push_back(dp);
    }
    double dt = 1.0 / frame_rate;
    double t = 0;

    vector<struct RectPoint> dp_points;

    struct RectPoint prev_trace;
    struct RectPoint curr_trace;

    struct Pixel trace_color;
    trace_color.r = 0x40;
    trace_color.g = 0x40;
    trace_color.b = 0x40;

    // #1e1e2e
    struct Pixel bg_color;
    bg_color.r = 0x1e;
    bg_color.g = 0x1e;
    bg_color.b = 0x2e;

    // #89b4fa
    struct Pixel pendulum_color;
    pendulum_color.r = 0x89;
    pendulum_color.g = 0xb4;
    pendulum_color.b = 0xfa;

    for (int frame = 0; frame < total_frames; frame++) {
        canvas.bg_fill(bg_color);

        for (int i = 0; i < num_pendulums; i++) {
            dp_arr[i]->update_state(t, dt);
            dp_points = dp_arr[i]->get_position(1);
            canvas.draw_line(origin, dp_points[0], pendulum_color);
            canvas.draw_line(dp_points[0], dp_points[1], pendulum_color);

            prev_trace = curr_trace;
            curr_trace = dp_points[1];
            if (frame != 0) { trace.draw_line(prev_trace, curr_trace, trace_color); }
        }
        t += dt;
        
        canvas += trace;
        
        string output = "frame_" + to_string(frame) + ".png";
        canvas.write_image(output);
    
    }
}

vector<double> get_perturbed(double di, double df, vector<double> perturbed_state, vector<double> initial_state) {
    vector<double> perturbed;
    for (int i = 0; i < perturbed_state.size(); i++) {
        perturbed.push_back(initial_state[i] + ((di * (perturbed_state[i] - initial_state[i])) / df));
    }
    return perturbed;
}

double temp_lyapunov(double di, double df) {
    return log(df / di);
}

double get_lyapunov_expo(double theta1, double theta2, double dt = 0.01, double run_time = 25, int sampling_freq = 20) {
    vector<double> lyapunov_exponents;

    DoublePendulum initial(theta1, theta2, 1, 1, 1, 1);
    DoublePendulum perturbed(theta1 + 0.001, theta2, 1, 1, 1, 1);
    double di = 0;
    for (int i = 0; i < initial.get_state().size(); i++) {
        di += (perturbed.get_state()[i] - initial.get_state()[i]) * (perturbed.get_state()[i] - initial.get_state()[i]);
    }
    di = sqrt(di);

    int counter = 0;

    for (double t = 0; t < run_time; t += dt) {
        counter++;
        initial.update_state(t, dt);
        perturbed.update_state(t, dt);

        if (!(counter % sampling_freq)) {
            double df = 0;
            for (int i = 0; i < initial.get_state().size(); i++) {
                df += (perturbed.get_state()[i] - initial.get_state()[i]) * (perturbed.get_state()[i] - initial.get_state()[i]);
            }
            df = sqrt(df);
            lyapunov_exponents.push_back(temp_lyapunov(di, df));
            vector<double> new_perturbed_state = get_perturbed(di, df, perturbed.get_state(), initial.get_state());
            perturbed.set_state(new_perturbed_state);
        }
    }

    double lyapunov_exponent = 0;
    for (int i = 0; i < lyapunov_exponents.size(); i++) {
        lyapunov_exponent += lyapunov_exponents[i];
    }
    return (lyapunov_exponent / run_time);
}

struct Pixel get_heatmap_color(double t) {
    const int NUM_COLORS = 4;
    static double color[NUM_COLORS][3] = {{0, 0, 1},  
                                          {0, 1, 0},  
                                          {1, 1, 0},  
                                          {1, 0, 0}};
    int idx1;      
    int idx2;      
    double frac_between = 0;
    
    if (t <= 0) {
        idx1 = idx2 = 0;            
    } else if (t >= 1) {
        idx1 = idx2 = NUM_COLORS - 1; 
    } else {
        t *= (NUM_COLORS - 1);      
        idx1 = t;                
        idx2 = idx1 + 1;                      
        frac_between = t - idx1;  
    }
    struct Pixel px;
        
    px.r = 255 * ((color[idx2][0] - color[idx1][0]) * frac_between + color[idx1][0]);
    px.g = 255 * ((color[idx2][1] - color[idx1][1]) * frac_between + color[idx1][1]);
    px.b = 255 * ((color[idx2][2] - color[idx1][2]) * frac_between + color[idx1][2]);

    return px;
}

void draw_hm_helper(uint8_t* img, int width, int start_height, int end_height, string thread_idx) {
    uint8_t *pix = img;

    double step = TAU / width;
    double theta1 = 0.0;
    double theta2 = step * start_height;

    for (int i = start_height; i < end_height; i++) {
        cout << "Thread " + thread_idx + ": " << ((i + 1 - start_height) / (double)(end_height - start_height)) * 100 << "%" << endl;
        theta2 += step;
        theta1 = 0;
        for (int j = 0; j < width; j++) {
            struct Pixel px;
            theta1 += step;
            //if ((((theta1 - M_PI) * (theta1 - M_PI)) + ((theta2 - M_PI) * (theta2 - M_PI))) < ((TAU / 3) * (TAU / 4))) {
                double norm_lya = get_lyapunov_expo(theta1, theta2) / 2.2;
                px = get_heatmap_color(norm_lya);
            //}
            pix = img + (((width * CHANNELS * j) + (i * CHANNELS)) - ((width * CHANNELS * j) + (i * CHANNELS)) % CHANNELS);
            *(pix + 0) = (px.r);
            *(pix + 1) = (px.g);
            *(pix + 2) = (px.b);
            pix += CHANNELS;
            //heatmap.draw_pixel(i, j, px);
        }
    }
}

void draw_heatmap_multithread(int width = 100, int height = 100) {
    if (height % NUM_THREADS != 0) {cout << "height % NUM_THREADS != 0" << endl; return;}
    int lines_per_thread = height / NUM_THREADS;

    int img_size = height * width * CHANNELS;
    uint8_t *img = new uint8_t[img_size];

    vector<thread*> thread_grp;

    for (int i = 0; i < NUM_THREADS; i++) {
        thread_grp.push_back(new thread(draw_hm_helper, ref(img), width, lines_per_thread * i, lines_per_thread * (i + 1), to_string(i)));
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        thread_grp[i]->join();
        delete thread_grp[i];
    }

    stbi_write_png("heatmap1.png", width, height, CHANNELS, img, width * CHANNELS);
    delete[] img;
}

void draw_heatmap_sing_thread(int width = 200, int height = 200) {
    const double TAU = 2 * M_PI;
    double step = TAU / (width);
    Graph heatmap(width, height, TAU, TAU);
    struct RectPoint init;

    for (int i = 0; i < height; i++) {
        //cout << "(" << init.x << ", " << init.y << ")" << endl;
        init.y += step;
        init.x = 0;
        for (int j = 0; j < width; j++) {
            init.x += step;
            double norm_lya = get_lyapunov_expo(init.x, init.y) / 2.2;
            struct Pixel px;
            px = get_heatmap_color(norm_lya);
            heatmap.draw_pixel(i, j, px);
        }
    }
    heatmap.write_image("heatmap.png");
}

vector<struct RectPoint> find_stable(const int& resolution, const double& threshold, const double& radius = M_PI_2) {
    double step = TAU / resolution;
    vector<struct RectPoint> stable;
    double min = 3;
    for (int i = 0; i < resolution; i++) {
        struct RectPoint init;
        init.x = i * step;
        for (int j = 0; j < resolution; j++) {
            init.y = j * step;
        
            if ((((init.x - M_PI) * (init.x - M_PI)) + ((init.y - M_PI) * (init.y - M_PI))) < (radius * radius)) {
                double lya = get_lyapunov_expo(init.x, init.y);
                if (lya < min) {
                    min = lya;
                    cout << min << endl;
                }
                if (lya < threshold) {
                    stable.push_back(init);
                    cout << "(" << setprecision(15) << init.x << ", " << init.y << ") -> " << lya << endl;
                }
            }
        }
    }
    return stable;
}

vector<struct RectPoint> stable_cluster(struct RectPoint& init, const int& resolution, const double& threshold, const double& radius = 0.1) {
    double step = 2 * radius / resolution;
    vector<struct RectPoint> stable;
    double min = 3;
    init.x -= radius;
    init.y -= radius;
    double inity = init.y;
    for (int i = 0; i < resolution; i++) {
        init.x += step;
        init.y = inity;
        for (int j = 0; j < resolution; j++) {
            init.y += step;
            double lya = get_lyapunov_expo(init.x, init.y);
            if (lya < min) {
                min = lya;
                cout << min << endl;
            }
            if (lya < threshold) {
                stable.push_back(init);
                cout << "(" << setprecision(15) << init.x << ", " << init.y << ") -> " << lya << endl;
            }
        }
    }
    return stable;
}

int main(int argc, char* argv[]) {
    //render(2.06780610206984, 2.46970906849588, 1, 1920, 1080, 3600);
    // 2.72533162698915, 4.43749962319558
    // 3.55576746631892, 1.84077694546277 < butterfly
    // 3.95460247116918, 3.01580622898317 < bee
    draw_heatmap_multithread(800, 800);
    //find_stable(2048, 0.2, (TAU / 4));
    // struct RectPoint init;
    // init.x = 3.55785;
    // init.y = 1.84569;
    // stable_cluster(init, 100, 0.3);

    // #1e1e2e
    // #89b4fa
    // #babbf1
}