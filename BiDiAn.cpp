/*
 * Bifurcation Diagram Animation v2.0
 * Lukas Kostal, 10.2.2023, ICL
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <format>
#include <filesystem>

using namespace std;

// define map to be animated
float map(float x, float r, float t) {
    x = exp(-t * pow(x, 2)) + r;
    return x;
}

// get order of number
int get_ord(int num) {
    int ord = 0;
    if (num > 0) {
        ord = floor(log10(num));
    }
    return ord;
}

int main() {

    // animation name
    string name = "Gauss_Map_3";

    // resolution
    int x_res = 1080; //4200
    int y_res = 720; //4200
    int inc = 2;

    // number of frames and fps
    int n_frm = 140;
    int fps = 20;

    // range of values of r
    float r_min = -1;
    float r_max = 0.8;

    // range of values of t
    float t_min = 3;
    float t_max = 13;

    // range of inital values of x
    float xi_min = 0;
    float xi_max = 1;

    // number of points per r and no of iterations
    int n_pts = 2000 ;
    int n_iter = 10;

    string frm_name;
    float prog;

    int k;
    float t;
    float r;
    float x;

    int n_vec = n_pts * x_res;
    vector<float> x_vec(n_vec, 0);
    vector<float> r_vec(n_vec, 0);

    float x_min;
    float x_max;

    int x_ind;
    int y_ind;

    int img_val;

    // print animation parameters
    cout << "animation:     " << name << endl;
    cout << "resolution:    " << x_res << "x" << y_res << endl;
    cout << "brightness:    " << inc << endl;
    cout << "no of frames:  " << n_frm << endl;
    cout << "framerate:     " << fps << endl;
    cout << endl;
    cout << "t:             " << "[" << t_min << ", " << t_max << "]" << endl;
    cout << "r:             " << "[" << r_min << ", " << r_max << "]" << endl;
    cout << "xi:            " << "[" << xi_min << ", " << xi_max << "]" << endl;
    cout << "points:        " << n_pts << endl;
    cout << "iterations:    " << n_iter << endl;
    cout << endl;

    // create or clear folder in home directory and set it as working directory
    string home_dir = (string) getenv("HOME");
    string folder_dir = home_dir + "/" + name;
    __fs::filesystem::remove_all(folder_dir);
    __fs::filesystem::create_directory(folder_dir);
    __fs::filesystem::current_path(folder_dir);

    for (int frm=0; frm < n_frm +1; frm++) {

        // generate alphanumeric frame name
        frm_name = "frm_";
        for (int i = 0; i < get_ord(n_frm) - get_ord(frm); i++) {
            frm_name += "0";
        }
        frm_name += to_string(frm);

        // update progress
        prog = (float)frm / (float)n_frm;
        cout << (float)(int)(prog * 1000) / 10 << "% \t [";

        for(int i=0; i < 50; i++) {
            if (i < ceil(prog * 50)){
                cout << "#";
            }
            else {
                cout << " ";
            }
        }
        cout << "] \t rendering " << frm_name << "\r" << flush;

        // initialise a 2D vector to hold image with all zeros
        vector<vector<uint8_t>> img_vec(x_res, vector<uint8_t> (y_res, 0));

        // find value of t for the current frame
        t = (t_max - t_min) * frm / n_frm + t_min;

        x_min = xi_max;
        x_max = xi_min;

        // loop over all values of r
        k = 0;
        for (int i=0; i < x_res; i++) {
            // calculate value of r
            r = (r_max - r_min) * i / x_res + r_min;

            for (int j=0; j < n_pts; j++) {
                // generate random inital value of x
                x = (float)rand() / (float)(RAND_MAX * (xi_max - xi_min)) + xi_min;

                // repeat for n iterations
                for (int l=0; l < n_iter; l++) {
                    x = map(x, r, t);
                }

                // write r and converged value of x to vector
                x_vec[k] = x;
                r_vec[k] = r;
                k ++;

                // find min and max values of converged x
                if (x < x_min) {
                    x_min = x;
                }
                if (x > x_max) {
                    x_max = x;
                }
            }
        }

        // index the image vector from values of x and r and increment pixel brightness
        for (int i=0; i < n_vec; i++) {
            x_ind = floor((r_vec[i] - r_min) * (x_res - 1) / (r_max - r_min));
            y_ind = floor((x_vec[i] - x_min) * (y_res - 1) / (x_max - x_min));
            if (img_vec[x_ind][y_ind] < 256 - inc) {
                    img_vec[x_ind][y_ind] += inc;
            }
        }

        // create the .ppm image file
        ofstream image;
        image.open(frm_name + ".ppm");

        // specify the properties of the image file to be generated
        image << "P3" << endl;
        image << x_res << " " << y_res << endl;
        image << "255" << endl;

        // write the image vector to the .ppm file
        for (int i = 0; i < y_res; i++) {
            for (int j = 0; j < x_res; j++) {
                img_val = +img_vec[j][y_res - i];
                image << img_val << " " << img_val << " " << img_val << endl;
            }
        }

        //close the image
        image.close();

        // call imagemagick to convert .ppm to .png and delete the .ppm
        system(("convert " + frm_name + ".ppm " + frm_name + ".png").c_str());
        std::__fs::filesystem::remove(frm_name + ".ppm");

    }

    // convert frames into .mp4 video with specified framerate
    system(("ffmpeg -r " + to_string(fps) + " -f image2 -s " + to_string(x_res) + "x" + to_string(y_res) + " -i frm_%0" \
            + to_string(get_ord(n_frm)+1) + "d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p" + name + ".mp4").c_str());
    return 0;
}
