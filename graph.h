/**
 * @file graph.cc
 * @author your name (you@domain.com)
 * @brief 
 */

#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <stdexcept>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace std;

const int CHANNELS = 3;

#ifndef STRUCT_GUARD
#define STRUCT_GUARD

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

#endif

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
    
    uint8_t *_img = nullptr;
    uint8_t *_pix = _img;
};

#endif