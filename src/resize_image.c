#include <math.h>
#include "image.h"
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable: 4996)
float nn_interpolate(image im, float x, float y, int c)
{
    // x : old coordinate col
    // y : old coordindate row
    // c : channel
    // im : old image
    int new_x_coord = round(x);
    int new_y_coord = round(y);
    
    float new_pixel = get_pixel(im,new_x_coord,new_y_coord,c);

    // TODO Fill in
    return new_pixel;
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    //float w_a, w_b;
    float w_a = (float)(im.w) / (w);
    float w_b = (float)0.5 * (w_a - 1);
    //float h_a, h_b;
    float h_a = (float)(im.h) / (h);
    float h_b = (float)0.5 * (h_a - 1);

    

    image resized = make_image(w, h, im.c);
    for (int row = 0; row < resized.h; row++) {
        for (int col = 0; col < resized.w; col++) {
            float old_x = w_a * col + w_b;
            float old_y = h_a * row + h_b;

            set_pixel(resized, col, row, 0, nn_interpolate(im,old_x,old_y,0));
            set_pixel(resized, col, row, 1, nn_interpolate(im, old_x, old_y, 1));
            set_pixel(resized, col, row, 2, nn_interpolate(im, old_x, old_y, 2));

        }
    }
    return resized;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    float new_pixel=0.;

    // TODO
    // 좌상단 모서리 밖일때
    //if (x < 0 && y < 0) {
    //    new_pixel = get_pixel(im, 0, 0, c);
    //}
    //// 우상단 모서리 밖일때
    //else if (y<0 && x>im.w-1) {
    //    new_pixel = get_pixel(im,im.w-1,0,c);

    //}
    //// 좌하단 모서리 밖일때
    //else if (x<0 && y>im.h - 1) {
    //    new_pixel = get_pixel(im, 0, im.h-1, c);

    //}
    //// 우하단 모서리 밖일때
    //else if (x>im.w-1 && y>im.h - 1) {
    //    new_pixel = get_pixel(im, im.w-1, im.h - 1, c);

    //}
    //// 맨위 edge 밖일때
    //else if (y < 0) {
    //    float v3_x = (float)floor(x);
    //    float v3_y = (float)0;
    //    float v4_x = (float)ceil(x);
    //    float v4_y = (float)0;
    //    float d1 = x - v3_x;
    //    float d2 = v4_x - x;
    //    new_pixel = get_pixel(im, v3_x, v3_y, c) * d2 + get_pixel(im, v4_x, v4_y, c) * d1;

    //}
    //// 맨아래 edge 밖일때
    //else if (y > im.h-1) {
    //    float v1_x = (float)floor(x);
    //    float v1_y = (float)im.h-1;
    //    float v2_x = (float)ceil(x);
    //    float v2_y = (float)im.h-1;
    //    float d1 = x - v1_x;
    //    float d2 = v2_x - x;
    //    new_pixel = get_pixel(im, v1_x, v1_y, c) * d2 + get_pixel(im, v2_x, v2_y, c) * d1;

    //}
    //// 맨왼쪽 edge 밖일때
    //else if (x < 0) {
    //    float v2_x = (float)0;
    //    float v2_y = (float)floor(y);
    //    float v4_x = (float)0;
    //    float v4_y = (float)ceil(y);
    //    float d3 = y-v2_y;
    //    float d4 = v4_y - y;
    //    new_pixel = get_pixel(im, v2_x, v2_y, c) * d4 + get_pixel(im, v4_x, v4_y, c) * d3;

    //}
    //// 맨오른쪽 edge 밖일때

    //else if (x >im.w-1) {
    //    float v1_x = (float)im.w-1;
    //    float v1_y = (float)floor(y);
    //    float v3_x = (float)im.w-1;
    //    float v3_y = (float)ceil(y);
    //    float d3 = y - v1_y;
    //    float d4 = v3_y - y;
    //    new_pixel = get_pixel(im, v1_x, v1_y, c) * d4 + get_pixel(im, v3_x, v3_y, c) * d3;

    //}
    // 일반형
    
    float v1_x = (float)floor(x);
    float v1_y = (float)floor(y);
    float v2_x = (float)ceil(x);
    float v2_y = (float)floor(y);
    float v3_x = (float)floor(x);
    float v3_y = (float)ceil(y);
    float v4_x = (float)ceil(x);
    float v4_y = (float)ceil(y);
    //int new_y_coord = round(y);
    //new_pixel=
    float d1 = x - v1_x;
    float d2 = v2_x - x;
    float d3 = y - v2_y;
    float d4 = v4_y - y;
    new_pixel = get_pixel(im, v1_x, v1_y, c) * d2 * d4 + get_pixel(im, v2_x, v2_y, c) * d1 * d4 + get_pixel(im, v3_x, v3_y, c) * d2 * d3 + get_pixel(im, v4_x, v4_y, c) * d1 * d3;
    
    return new_pixel;
}

image bilinear_resize(image im, int w, int h)
{
    float w_a = (float)(im.w) / (w);
    float w_b = (float)0.5 * (w_a - 1);
    //float h_a, h_b;
    float h_a = (float)(im.h) / (h);
    float h_b = (float)0.5 * (h_a - 1);



    image resized = make_image(w, h, im.c);
    for (int row = 0; row < resized.h; row++) {
        for (int col = 0; col < resized.w; col++) {
            float old_x = w_a * col + w_b;
            float old_y = h_a * row + h_b;

            set_pixel(resized, col, row, 0, bilinear_interpolate(im, old_x, old_y, 0));
            set_pixel(resized, col, row, 1, bilinear_interpolate(im, old_x, old_y, 1));
            set_pixel(resized, col, row, 2, bilinear_interpolate(im, old_x, old_y, 2));

        }
    }
    return resized;
}

