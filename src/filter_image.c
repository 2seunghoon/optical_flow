#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable: 4996)
void l1_normalize(image im)
{
    // TODO
    for (int channel = 0; channel < im.c; channel++) {
        float pixel_sum = 0.;

        for (int row = 0; row < im.h; row++) {
            for (int col = 0; col < im.w; col++) {
                pixel_sum += get_pixel(im, col, row, channel);
            }
        }
        for (int row = 0; row < im.h; row++) {
            for (int col = 0; col < im.w; col++) {

                set_pixel(im, col, row, channel, get_pixel(im, col, row, channel) / pixel_sum);
            }
        }
    }
}

image make_box_filter(int w)
{
    image box_filter;
    // TODO
    box_filter = make_image(w, w, 1);
    for (int i = 0; i < w * w; ++i) {
        box_filter.data[i] = 1.0;
    }
    l1_normalize(box_filter);

    return box_filter;
}


image convolve_image(image im, image filter, int preserve)
{
    // TODO
    image new_image = make_image(im.w,im.h,im.c);
    if (preserve) {
        for (int c = 0; c < im.c; c++) {
            for (int row = 0; row < im.h; row++) {
                for (int col = 0; col < im.w; col++) {
                    float weight_sum = 0.;
                    int x0 = col - filter.w / 2;
                    int y0 = row- filter.h / 2;

                    for (int filter_row = 0; filter_row < filter.h; filter_row++) {
                        for (int filter_col = 0; filter_col < filter.w; filter_col++) {
                            weight_sum += get_pixel(im, x0 + filter_row, y0 + filter_col, c) * get_pixel(filter, filter_row, filter_col, 0);
                        }
                    }
                    set_pixel(new_image, col, row, c, weight_sum);

                }
            }
        }
    }
    else {
        new_image = make_image(im.w, im.h, 1);
        for (int row = 0; row < im.h; row++) {
            for (int col = 0; col < im.w; col++) {
                float weight_sum = 0.;
                int x0 = col - filter.w / 2;
                int y0 = row - filter.h / 2;

                for (int filter_row = 0; filter_row < filter.h; filter_row++) {
                    for (int filter_col = 0; filter_col < filter.w; filter_col++) {
                        if (im.c == 3) {
                            weight_sum += get_pixel(im, x0 + filter_row, y0 + filter_col, 0) * get_pixel(filter, filter_row, filter_col, 0);
                            weight_sum += get_pixel(im, x0 + filter_row, y0 + filter_col, 1) * get_pixel(filter, filter_row, filter_col, 0);
                            weight_sum += get_pixel(im, x0 + filter_row, y0 + filter_col, 2) * get_pixel(filter, filter_row, filter_col, 0);
                        }
                        else {
                            weight_sum += get_pixel(im, x0 + filter_row, y0 + filter_col, 0) * get_pixel(filter, filter_row, filter_col, 0);

                        }
                    }
                }
                set_pixel(new_image, col, row, 0, weight_sum);

            }
        }
        

    }
    return new_image;
}

image make_highpass_filter()
{
    // TODO
    image highpass = make_image(3, 3, 1);

    set_pixel(highpass, 0, 0, 0, 0.);
    set_pixel(highpass, 1, 0, 0, -1.);
    set_pixel(highpass, 2, 0, 0, 0.);
    set_pixel(highpass, 0, 1, 0, -1.);
    set_pixel(highpass, 1, 1, 0, 4.);
    set_pixel(highpass, 2, 1, 0, -1.);
    set_pixel(highpass, 0, 2, 0, 0.);
    set_pixel(highpass, 1, 2, 0, -1.);
    set_pixel(highpass, 2, 2, 0, 0.);


    /*highpass.data[0] = 0.;
    highpass.data[1] = -1.;
    highpass.data[2] = 0.;
    highpass.data[3] = -1.;
    highpass.data[4] = 4.;
    highpass.data[5] = -1.;
    highpass.data[6] = 0.;
    highpass.data[7] = -1.;
    highpass.data[8] = 0.;*/

    return highpass;
}

image make_sharpen_filter()
{
    // TODO
    
    image sharpen = make_image(3, 3, 1);
    sharpen.data[0] = 0.;
    sharpen.data[1] = -1.;
    sharpen.data[2] = 0.;
    sharpen.data[3] = -1.;
    sharpen.data[4] = 5.;
    sharpen.data[5] = -1.;
    sharpen.data[6] = 0.;
    sharpen.data[7] = -1.;
    sharpen.data[8] = 0.;
    return sharpen;
}

image make_emboss_filter()
{
    // TODO
    image emboss = make_image(3, 3, 1);
    emboss.data[0] = -2.;
    emboss.data[1] = -1.;
    emboss.data[2] = 0.;
    emboss.data[3] = -1.;
    emboss.data[4] = 1.;
    emboss.data[5] = 1.;
    emboss.data[6] = 0.;
    emboss.data[7] = 1.;
    emboss.data[8] = 2.;
    return emboss;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    int size_ = sigma * 6.;
    int size;
    if (size_ % 2 == 0) {
        size = size_ + 1;
    }
    else {
        size = size_ + 2;

    }
    int x0 = size / 2;
    int y0 = x0;
    image gaussian = make_image(size, size, 1);
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            float values = 1 / (TWOPI*sigma*sigma)*exp(-((row-x0)*(row-x0)+(col-y0)*(col-y0))/(2*sigma*sigma));
            set_pixel(gaussian, col, row, 0, values);
        }
    }
    l1_normalize(gaussian);
    return gaussian;
}

image add_image(image a, image b)
{
    // TODO
    image new_image = make_image(a.w, a.h, a.c);
    for (int i = 0; i < a.w * a.h * a.c; i++) {
        new_image.data[i] = a.data[i] + b.data[i];
    }
    return new_image;
}

image sub_image(image a, image b)
{
    // TODO
    image new_image = make_image(a.w, a.h, a.c);
    for (int i = 0; i < a.w * a.h * a.c; i++) {
        new_image.data[i] = a.data[i] - b.data[i];
    }
    return new_image;
}

image make_gx_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);
    filter.data[0] = -1.;
    filter.data[1] = 0.;
    filter.data[2] = 1.;
    filter.data[3] = -2.;
    filter.data[4] = 0.;
    filter.data[5] = 2.;
    filter.data[6] = -1.;
    filter.data[7] = 0.;
    filter.data[8] = 1.;

    return filter;
}

image make_gy_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);
    filter.data[0] = -1.;
    filter.data[1] = -2.;
    filter.data[2] = -1.;
    filter.data[3] = 0.;
    filter.data[4] = 0.;
    filter.data[5] = 0.;
    filter.data[6] = 1.;
    filter.data[7] = 2.;
    filter.data[8] = 1.;
    return filter;
}

void feature_normalize(image im)
{
    // TODO
    // minmax scaler
    float max_ = im.data[0];
    float min_ = im.data[0];
    for (int i = 0; i < im.w * im.h * im.c; i++) {
        max_ = fmaxf(max_, im.data[i]);
        min_ = fminf(min_, im.data[i]);

    }
    float range = max_ - min_;
    if (range == 0.) {
        for (int i = 0; i < im.w * im.h * im.c; i++) {
            im.data[i] = 0.;
        }
    }
    else {
        for (int i = 0; i < im.w * im.h * im.c; i++) {
            im.data[i] = (im.data[i]-min_)/(range);
        }
    }

}

image *sobel_image(image im)
{
    // TODO
    // 0 for mag, 1 for direction
    image* res = calloc(2, sizeof(image));
    //res[0].w = im.w;
    //res[0].h = im.h;
    //res[0].c = 0; 
    res[0] = make_image(im.w, im.h, 1);
    res[1] = make_image(im.w, im.h, 1);

    //res[1].w = im.w;
    //res[1].h = im.h;
    //res[1].c = 0;

    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image gx = convolve_image(im, gx_filter, 0);
    image gy = convolve_image(im, gy_filter, 0);

    for (int i = 0; i < im.w * im.h; i++) {
        res[0].data[i] = sqrtf(gx.data[i] * gx.data[i] + gy.data[i] * gy.data[i]);
        res[1].data[i] = atan2f(gy.data[i], gx.data[i]);
    }
    return res;
}

image colorize_sobel(image im)
{
    // TODO
    image* res = sobel_image(im);
    image new_image = make_image(im.w, im.h, im.c);
    image mag = res[0];
    image dir = res[1];
    feature_normalize(mag);
    feature_normalize(dir);

    for (int row = 0; row < im.h; row++) {
        for (int col = 0; col < im.w; col++) {
            set_pixel(new_image, col, row, 0, get_pixel(dir, col, row, 0));
            set_pixel(new_image, col, row, 1, get_pixel(mag, col, row, 0));
            set_pixel(new_image, col, row, 2, get_pixel(mag, col, row, 0));

        }
    }
    hsv_to_rgb(new_image);
    return convolve_image(new_image,make_gaussian_filter(1),1);
}
