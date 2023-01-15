#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable: 4996)
float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    //return im.data[c][x][y]
    if (im.w <= x) {
        x = im.w - 1;
    }
    else if (x < 0) {
        x = 0;
    }
    if (im.h <= y) {
        y = im.h- 1;
    }
    else if (y < 0) {
        y = 0;
    }
    float pixel_value = im.data[c * im.w * im.h + im.w * y + x];
    return pixel_value;
    //return 0;
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if (im.w <= x) {
        x = im.w - 1;
    }
    else if (x < 0) {
        x = 0;
    }
    if (im.h <= y) {
        y = im.h - 1;
    }
    else if (y < 0) {
        y = 0;
    }
    im.data[c * im.w * im.h + im.w * y + x] = v;
    
    // TODO Fill this in
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    
    memcpy(copy.data, im.data,im.w*im.h*im.c*sizeof(float));
    //copy.data = im.data;
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    int index = 0;
    for (int row = 0; row < im.h; row++) {
        for (int col = 0; col < im.w; col++) {
            float r_pixel = get_pixel(im, col, row, 0);
            float g_pixel = get_pixel(im, col, row, 1);
            float b_pixel = get_pixel(im, col, row, 2);

            float gray_pixel = 0.299 * r_pixel + 0.587 * g_pixel + 0.114 * b_pixel;
            gray.data[index++] = gray_pixel;

        }
    }
    
    
    return gray;
}

void shift_image(image im, int c, float v)
{
    for (int row = 0; row < im.h ; row++) {
        for (int col = 0; col < im.w; col++) {
            float pixel = get_pixel(im, col, row, c);
            pixel += v;
            /*if (pixel >= 1.0) {
                pixel = 1.0;
            }
            else if (pixel < 0) {
                pixel = 0;
            }*/
            set_pixel(im, col, row, c, pixel);
        }
    }
    // TODO Fill this in
}

void clamp_image(image im)
{
    for (int index = 0; index < im.w * im.c * im.h; index++) {
        if (im.data[index] >= 1.0) {
            im.data[index] = 1.0;
        }
        else if (im.data[index] < 0) {
            im.data[index] = 0;
        }
    }
    // TODO Fill this in
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    for (int row = 0; row < im.h; row++) {
        for (int col = 0; col < im.w; col++) {
            float V = three_way_max(get_pixel(im, col, row, 0), get_pixel(im, col, row, 1), get_pixel(im, col, row, 2));
            float S = 0.;
            float m = three_way_min(get_pixel(im, col, row, 0), get_pixel(im, col, row, 1), get_pixel(im, col, row, 2));
            float C = V - m;
            if (V == 0) {
                S = 0;
            }
            else {
                S = C / V;
            }

            float H = 0.;
            float H_prime = 0.;
            if (C == 0) {
                H = 0;
            }
            else {
                float R_pixel = get_pixel(im, col, row, 0);
                float G_pixel = get_pixel(im, col, row, 1);
                float B_pixel = get_pixel(im, col, row, 2);

                if (V == R_pixel) {
                    H_prime = (G_pixel - B_pixel) / C;
                }
                else if (V == G_pixel) {
                    H_prime = (B_pixel - R_pixel) / C + 2;
                }
                else if (V == B_pixel) {
                    H_prime = (R_pixel - G_pixel) / C +4;

                }

                if (H_prime < 0) {
                    H = H_prime / 6 + 1;
                }
                else {
                    H = H_prime / 6;
                }
            }
            set_pixel(im, col, row, 0, H);
            set_pixel(im, col, row, 1, S);
            set_pixel(im, col, row, 2, V);


        }
    }
    // TODO Fill this in
}

void hsv_to_rgb(image im)
{
    for (int row = 0; row < im.h; row++) {
        for (int col = 0; col < im.w; col++) {
            float H = get_pixel(im, col, row, 0);
            float S = get_pixel(im, col, row, 1);
            float V = get_pixel(im, col, row, 2);
            H *= 6.; // for H=[0,6)
            
            float M = V;
            //float m = M*(1 - S);
            float C = V * S;
            float X = (C) * (1 - fabs(fmod(H,2) - 1));
            float m = V - C;
            float r, g, b;
            /*if (H == 0. ||V==0||S==0) {
                r = 0;
                g = 0;
                b = 0;
            }*/
            if (H >= 0 && H < 1) {
                //r = C, g = X, b = 0;
                r = C+m, g = X+m, b = m;

            }
            else if (H >= 1 && H < 2) {
                //r = X, g = C, b = 0;
                r = X+m, g = C+m, b =m;

            }
            else if (H >= 2 && H < 3) {
                r = m, g = C+m, b = X+m;
            }
            else if (H >= 3 && H < 4) {
                r = m, g = X+m, b = C+m;
            }
            else if (H >= 4 && H < 5) {
                r = X+m, g = m, b = C+m;
            }
            else if (H >= 5 && H < 6) {
                r = C+m, g =m, b = X+m;
            }
            else {
                r = m;
                g = m;
                b = m;
            }

            set_pixel(im, col, row, 0, r);
            set_pixel(im, col, row, 1, g);
            set_pixel(im, col, row, 2, b);
        }
    }
    // TODO Fill this in
}
