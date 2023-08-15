#include <vector>
#include <cmath>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model *model = nullptr;
const int width = 800;
const int height = 800;

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color)
{
    bool steep = false;
    if (std::abs(x0 - x1) < std::abs(y0 - y1))
    { // if the line is steep, we transpose the image
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0 > x1)
    { // make it left−to−right
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int dx = x1 - x0;
    int dy = y1 - y0;
    int derror2 = std::abs(dy) * 2;
    int error2 = 0;
    int y = y0;
    for (int x = x0; x <= x1; x++)
    {
        if (steep)
        {
            image.set(y, x, color); // if transposed, de−transpose
        }
        else
        {
            image.set(x, y, color);
        }
        error2 += derror2;
        if (error2 > dx)
        {
            y += (y1 > y0 ? 1 : -1);
            error2 -= dx * 2;
        }
    }
}

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color)
{
    bool steep = false;
    if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y))
    {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    if (p0.x > p1.x)
    {
        std::swap(p0, p1);
    }

    for (int x = p0.x; x <= p1.x; x++)
    {
        float t = (x - p0.x) / (float)(p1.x - p0.x);
        int y = p0.y * (1. - t) + p1.y * t;
        if (steep)
        {
            image.set(y, x, color);
        }
        else
        {
            image.set(x, y, color);
        }
    }
}

// Vec3f barycentric(Vec2i *pts, Vec2i P)
// {
//     Vec3f u = Vec3f(pts[2][0] - pts[0][0], pts[1][0] - pts[0][0], pts[0][0] - P[0]) ^ Vec3f(pts[2][1] - pts[0][1], pts[1][1] - pts[0][1], pts[0][1] - P[1]);
//     /* `pts` and `P` has integer value as coordinates
//        so `abs(u[2])` < 1 means `u[2]` is 0, that means
//        triangle is degenerate, in this case return something with negative coordinates */
//     if (std::abs(u.z) < 1)
//         return Vec3f(-1, 1, 1);
//     return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
// }

Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
{
    Vec3f s[2];
    for (int i = 2; i--;)
    {
        s[i][0] = C[i] - A[i];
        s[i][1] = B[i] - A[i];
        s[i][2] = A[i] - P[i];
    }
    Vec3f u = cross(s[0], s[1]);
    if (std::abs(u[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
    return Vec3f(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

// void triangle(Vec2i *pts, TGAImage &image, TGAColor color)
// {
//     Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
//     Vec2i bboxmax(0, 0);
//     Vec2i clamp(image.get_width() - 1, image.get_height() - 1);
//     for (int i = 0; i < 3; i++)
//     {
//         bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
//         bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));

//         bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
//         bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
//     }
//     Vec2i P;
//     for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++)
//     {
//         for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++)
//         {
//             Vec3f bc_screen = barycentric(pts, P);
//             if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0)
//                 continue;
//             image.set(P.x, P.y, color);
//         }
//     }
// }

void triangle(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color)
{
    Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            bboxmin[j] = std::max(0.f, std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++)
    {
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++)
        {
            Vec3f bc_screen = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0)
                continue;
            P.z = 0;
            for (int i = 0; i < 3; i++)
                P.z += pts[i][2] * bc_screen[i];
            if (zbuffer[int(P.x + P.y * width)] < P.z)
            {
                zbuffer[int(P.x + P.y * width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}

Vec3f world2screen(Vec3f v)
{
    return Vec3f(int((v.x + 1.) * width / 2. + .5), int((v.y + 1.) * height / 2. + .5), v.z);
}

int main(int argc, char **argv)
{
    if (2 == argc)
    {
        model = new Model(argv[1]);
    }
    else
    {
        model = new Model("obj/african_head.obj");
    }
    // TGAImage image(width, height, TGAImage::RGB);
    // // line(13, 20, 80, 40, image, white);
    // for (int i = 0; i < model->nfaces(); i++)
    // {
    //     std::vector<int> face = model->face(i);
    //     for (int j = 0; j < 3; j++)
    //     {
    //         Vec3f v0 = model->vert(face[j]);
    //         Vec3f v1 = model->vert(face[(j + 1) % 3]);
    //         int x0 = (v0.x + 1.) * width / 2.;
    //         int y0 = (v0.y + 1.) * height / 2.;
    //         int x1 = (v1.x + 1.) * width / 2.;
    //         int y1 = (v1.y + 1.) * height / 2.;
    //         line(x0, y0, x1, y1, image, white);
    //     }
    // }

    // Vec2i t0[3] = {Vec2i(10, 70), Vec2i(50, 160), Vec2i(70, 80)};
    // Vec2i t1[3] = {Vec2i(180, 50), Vec2i(150, 1), Vec2i(70, 180)};
    // Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)};

    // triangle(t0[0], t0[1], t0[2], image, red);
    // triangle(t1[0], t1[1], t1[2], image, white);
    // triangle(t2[0], t2[1], t2[2], image, green);

    TGAImage frame(800, 800, TGAImage::RGB);
    // Vec2i pts[3] = {Vec2i(10, 10), Vec2i(100, 30), Vec2i(190, 160)};
    // triangle(pts, frame, TGAColor(255, 0, 0));

    Vec3f light_dir(0, 0, -1); // define light_dir

    float *zbuffer = new float[width * height];
    for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max())
        ;

    for (int i = 0; i < model->nfaces(); i++)
    {
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        Vec3f world_coords[3];
        for (int j = 0; j < 3; j++)
        {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
            world_coords[j] = v;
        }
        Vec3f n = cross((world_coords[2] - world_coords[0]), (world_coords[1] - world_coords[0]));
        n.normalize();
        float intensity = n * light_dir;
        Vec3f pts[3];
        for (int i = 0; i < 3; i++)
            pts[i] = world2screen(model->vert(face[i]));
        if (intensity > 0)
        {
            triangle(pts, zbuffer, frame, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
        }
    }

    frame.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    frame.write_tga_file("output.tga");
    delete model;
    return 0;
}