/**
 * @file graph.h
 * @author Garrett
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
#include "structs.h"


using namespace std;

const int CHANNELS = 3;

class Graph {
public:
    Graph();
    Graph(const int& img_width_res,  const int& img_height_res, 
          const double& graph_width, const double& graph_height);
    ~Graph();

    int get_img_width() const;
    int get_img_height() const;
    double get_graph_width() const;
    double get_graph_height() const;
    int get_img_size() const;
    uint8_t* get_img() const;

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

    //void write_image(const string& file_name);

    Graph& operator-=(Graph& graph);
    Graph& operator+=(Graph& graph);
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