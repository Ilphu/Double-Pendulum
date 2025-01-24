/**
 * @file main.cc
 * @author Garrett Rhoads
 * @brief A study in stable double pendulum orbits
 * @date 2024-12-25
 */

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <stdexcept>
//#include <iomapip>
//#include <bits/stdc++.h> 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "graph.h"
#include "double_pendulum.h"
#include "structs.h"

#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))


using namespace std;

const int NUM_THREADS = 4;
//const int CHANNELS = 3;
const double TAU = 2 * M_PI;

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
        
        string output = "frames/frame_" + to_string(frame) + ".png";
        stbi_write_png(output.c_str(), canvas.get_img_width(), 
                       canvas.get_img_height(), CHANNELS, canvas.get_img(), 
                       canvas.get_img_width() * CHANNELS);

        for (unsigned int i = 0; i < dp_arr.size(); i++) {
            delete dp_arr[i];
        }
    }
}

vector<double> get_perturbed(double di, double df, vector<double> perturbed_state, 
                                                     vector<double> initial_state) {
    vector<double> perturbed;
    for (unsigned int i = 0; i < perturbed_state.size(); i++) {
        perturbed.push_back(initial_state[i] + ((di * (perturbed_state[i] - initial_state[i])) / df));
    }
    return perturbed;
}

double temp_lyapunov(double di, double df) {
    return log(df / di);
}

double get_lyapunov_expo(double theta1, double theta2, double dt = 0.01, 
                         double run_time = 25, int sampling_freq = 20) {
    vector<double> lyapunov_exponents;

    DoublePendulum initial(theta1, theta2, 1, 1, 1, 1);
    DoublePendulum perturbed(theta1 + 0.001, theta2, 1, 1, 1, 1);
    double di = 0;
    for (unsigned int i = 0; i < initial.get_state().size(); i++) {
        di += (perturbed.get_state()[i] - initial.get_state()[i]) * 
              (perturbed.get_state()[i] - initial.get_state()[i]);
    }
    di = sqrt(di);

    int counter = 0;

    for (double t = 0; t < run_time; t += dt) {
        counter++;
        initial.update_state(t, dt);
        perturbed.update_state(t, dt);

        if (!(counter % sampling_freq)) {
            double df = 0;
            for (unsigned int i = 0; i < initial.get_state().size(); i++) {
                df += (perturbed.get_state()[i] - initial.get_state()[i]) * 
                      (perturbed.get_state()[i] - initial.get_state()[i]);
            }
            df = sqrt(df);
            lyapunov_exponents.push_back(temp_lyapunov(di, df));
            vector<double> new_perturbed_state = get_perturbed(di, df, perturbed.get_state(), initial.get_state());
            perturbed.set_state(new_perturbed_state);
        }
    }

    double lyapunov_exponent = 0;
    for (unsigned int i = 0; i < lyapunov_exponents.size(); i++) {
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

void draw_hm_helper(uint8_t* img, int width, int height, int start_height, int end_height, string thread_idx) {
    uint8_t *pix = img;

    double step = TAU / width;
    double theta1 = step * start_height;
    double theta2 = 0.0;

    for (int i = start_height; i < end_height; i++) {
        cout << "Thread " + thread_idx + ": " << ((i + 1 - start_height) / 
                    (double)(end_height - start_height)) * 100 << "%" << endl;
        theta1 += step;
        theta2 = 0;
        for (int j = 0; j < width; j++) {
            struct Pixel px;
            theta2 += step;
            //if ((((theta1 - M_PI) * (theta1 - M_PI)) + ((theta2 - M_PI) * (theta2 - M_PI))) < ((TAU / 3) * (TAU / 4))) {
                double norm_lya = get_lyapunov_expo(theta1, theta2) / 2.2;
                px = get_heatmap_color(norm_lya);
            //}
            pix = img + (((width * CHANNELS * i) + (j * CHANNELS)) - 
                         ((width * CHANNELS * i) + (j * CHANNELS)) % CHANNELS);
            *(pix + 0) = (px.r);
            *(pix + 1) = (px.g);
            *(pix + 2) = (px.b);
            pix = img + (((width * CHANNELS * (height - i - 1)) + ((width - j - 1) * CHANNELS)) - 
                         ((width * CHANNELS * (height - i - 1)) + ((width - j - 1) * CHANNELS)) % CHANNELS);
            *(pix + 0) = (px.r);
            *(pix + 1) = (px.g);
            *(pix + 2) = (px.b);
            //heatmap.draw_pixel(i, j, px);
        }
    }
}

void draw_heatmap_multithread(int width = 100, int height = 100) {
    if (height % NUM_THREADS != 0) {cout << "height % NUM_THREADS != 0" << endl; return;}
    int lines_per_thread = height / (NUM_THREADS * 2);

    int img_size = height * width * CHANNELS;
    uint8_t *img = new uint8_t[img_size];

    vector<thread*> thread_grp;

    for (int i = 0; i < NUM_THREADS; i++) {
        thread_grp.push_back(new thread(draw_hm_helper, ref(img), width, height, 
                lines_per_thread * i, lines_per_thread * (i + 1), to_string(i)));
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
    string output = "heatmap.png";
    stbi_write_png(output.c_str(), heatmap.get_img_width(), 
                   heatmap.get_img_height(), CHANNELS, heatmap.get_img(), 
                   heatmap.get_img_width() * CHANNELS);
    //heatmap.write_image("heatmap.png");
}

vector<struct RectPoint> find_stable(const int& resolution, const double& threshold, 
                                                        const double& radius = M_PI_2) {
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
                    //cout << "(" << setprecision(15) << init.x << ", " << init.y << ") -> " << lya << endl;
                }
            }
        }
    }
    return stable;
}

vector<struct RectPoint> stable_cluster(struct RectPoint& init, const int& resolution, 
                                        const double& threshold, const double& radius = 0.1) {
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
                //cout << "(" << setprecision(15) << init.x << ", " << init.y << ") -> " << lya << endl;
            }
        }
    }
    return stable;
}

int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
        if (argv[i] == "-r") {
            render(2.06780610206984, 2.46970906849588, 1, 1920, 1080, 300);
        } else if (argv[i] == "-hm") {
            draw_heatmap_multithread(64, 64);
        }
    }
    //render(2.06780610206984, 2.46970906849588, 1, 1920, 1080, 300);
    // 2.72533162698915, 4.43749962319558
    // 3.55576746631892, 1.84077694546277 < butterfly
    // 3.95460247116918, 3.01580622898317 < bee
    //draw_heatmap_multithread(64, 64);
    //find_stable(2048, 0.2, (TAU / 4));
    // struct RectPoint init;
    // init.x = 3.55785;
    // init.y = 1.84569;
    // stable_cluster(init, 100, 0.3);
}