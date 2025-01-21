/**
 * @file graph.cc
 * @author your name (you@domain.com)
 * @brief 
 */

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <stdexcept>
#include "graph.h"

using namespace std;

Graph::Graph() {}

Graph::Graph(const int& img_width_res, const int& img_height_res, 
             const double& graph_width, const double& graph_height) {
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