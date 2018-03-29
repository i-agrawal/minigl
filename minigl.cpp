/**
 * Name: Ishan Agrawal
 * SID:  861229134
 * user: iagra001
 *
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

#define POS 0
#define COL 1
typedef vector<vec4> vertex;

#define A 0
#define B 1
#define C 2
typedef vector<vertex> triangle;

/**
 * Global data types for usage
 */
vec4 color;
MGLpoly_mode pmode;
vector<vertex> points;
vector<triangle> triangles;
MGLmatrix_mode mmode;
vector<mat4> stacks[2] = {{mat4()},{mat4()}};
MGLfloat* zbuffer = 0;


/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
vec3 px(vec4 pt, MGLsize w, MGLsize h) {
    int x = w*(pt[0]+1)/2;
    int y = h*(pt[1]+1)/2;
    return vec3(x-0.5,y-0.5,pt[2]);
}

float area(const vec3 &a, const vec3 &b, const vec3 &c)
{
    return a[0]*(b[1]-c[1])+a[1]*(c[0]-b[0])+(b[0]*c[1]-b[1]*c[0]);
}

bool clipped(vec3 v)
{
    return (v[0] < -1 || v[0] > 1 || v[1] < -1 || v[1] > 1 || v[2] < -1 || v[2] > 1);
}

void Rasterize_Triangle(const triangle &t, MGLsize w, MGLsize h, MGLpixel *data)
{
    vec3 a = px(t[A][POS],w,h);
    vec3 b = px(t[B][POS],w,h);
    vec3 c = px(t[C][POS],w,h);

    int minx = min(a[0],min(b[0],c[0]));
    int miny = min(a[1],min(b[1],c[1]));
    int maxx = max(a[0],max(b[0],c[0]));
    int maxy = max(a[1],max(b[1],c[1]));

    float abc = area(a,b,c);

    int i, j;
    for (i = max(minx,0); i <= min(maxx,(int)w-1); i++) {
        for (j = max(miny,0); j <= min(maxy,(int)h-1); j++) {
            vec3 p(i,j,0);
            float pbc = area(p,b,c) / abc;
            float apc = area(a,p,c) / abc;
            float abp = 1 - apc - pbc;

            if (pbc >= 0 && apc >= 0 && abp >= 0) {
                double x = 2.0*i/w - 1.0;
                double y = 2.0*j/h - 1.0;
                double z = pbc*t[A][POS][2] + apc*t[B][POS][2] + abp*t[C][POS][2];
                if (!clipped(vec3(x,y,z)) && z < zbuffer[i+j*w]) {
                    zbuffer[i+j*w] = z;
                    vec4 pxcolor = t[A][COL]*pbc + t[B][COL]*apc + t[C][COL]*abp;
                    data[i+j*w] = Make_Pixel(pxcolor[0]*255,
                                             pxcolor[1]*255,
                                             pxcolor[2]*255);
                }
            }
        }
    }
}

void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
    if (zbuffer != 0) free(zbuffer);
    zbuffer = (MGLfloat*)malloc(width*height*sizeof(MGLfloat));

    size_t i;
    for (i = 0; i < width*height; i++) {
        data[i] = Make_Pixel(0,0,0);
        zbuffer[i] = 2;
    }

    for (i = 0; i < triangles.size(); i++)
        Rasterize_Triangle(triangles[i], width, height, data);

    triangles.clear();
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
    pmode = mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
    size_t i;
    switch(pmode) {
        case MGL_TRIANGLES:
            for (i = 0; i < points.size(); i+=3)
                triangles.push_back({points[i], points[i+1], points[i+2]});
            break;
        case MGL_QUADS:
            for (i = 0; i < points.size(); i+=4) {
                triangles.push_back({points[i], points[i+1], points[i+2]});
                triangles.push_back({points[i], points[i+3], points[i+2]});
            }
            break;
        default:
            MGL_ERROR("GL_INVALID_ENUM");
            break;
    }

    points.clear();
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
    mglVertex3(x,y,0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
    vec4 pt(x,y,z,1);
    pt = stacks[MGL_MODELVIEW].back() * pt;
    pt = stacks[MGL_PROJECTION].back() * pt;
    pt /= pt[3];
    points.push_back({pt,color});
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    mmode = mode;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
    mat4 cp;
    for (int i = 0; i < 16; i++)
        cp.values[i] = stacks[mmode].back().values[i];

    stacks[mmode].push_back(cp);
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
    stacks[mmode].pop_back();
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
    for (int i = 0; i < 16; i++)
        stacks[mmode].back().values[i] = (i%5 == 0);
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
    for (int i = 0; i < 16; i++)
        stacks[mmode].back().values[i] = matrix[i];
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
    mat4 matx;
    for (int i = 0; i < 16; i++)
        matx.values[i] = matrix[i];

    stacks[mmode].back() = stacks[mmode].back() * matx;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
    mat4 trans = {{0}};

    trans(0,0) = 1;
    trans(0,3) = x;
    trans(1,1) = 1;
    trans(1,3) = y;
    trans(2,2) = 1;
    trans(2,3) = z;
    trans(3,3) = 1;

    mglMultMatrix(trans.values);
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
    float c = cos(angle * M_PI / 180);
    float s = sin(angle * M_PI / 180);
    vec3 orig(x,y,z);
    orig = orig.normalized();
    x = orig[0]; y = orig[1]; z = orig[2];


    mat4 rotate = {{0}};
    rotate(0,0) = x*x*(1-c) + c;
    rotate(0,1) = x*y*(1-c) - z*s;
    rotate(0,2) = x*z*(1-c) + y*s;
    rotate(1,0) = y*x*(1-c) + z*s;
    rotate(1,1) = y*y*(1-c) + c;
    rotate(1,2) = y*z*(1-c) - x*s;
    rotate(2,0) = z*x*(1-c) - y*s;
    rotate(2,1) = z*y*(1-c) + x*s;
    rotate(2,2) = z*z*(1-c) + c;
    rotate(3,3) = 1;

    mglMultMatrix(rotate.values);
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
    mat4 scale = {{0}};
    scale(0,0) = x;
    scale(1,1) = y;
    scale(2,2) = z;
    scale(3,3) = 1;

    mglMultMatrix(scale.values);
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
    mat4 frust = {{0}};
    frust(0,0) = 2*near/(right-left);
    frust(0,2) = (right+left)/(right-left);
    frust(1,1) = 2*near/(top-bottom);
    frust(1,2) = (top+bottom)/(top-bottom);
    frust(2,2) = (far+near)/(near-far);
    frust(2,3) = (2*far*near)/(near-far);
    frust(3,2) = -1;

    mglMultMatrix(frust.values);
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
    mat4 ortho = {{0}};
    ortho(0,0) = 2/(right-left);
    ortho(0,3) = (right+left)/(left-right);
    ortho(1,1) = 2/(top-bottom);
    ortho(1,3) = (top+bottom)/(bottom-top);
    ortho(2,2) = -2/(far-near);
    ortho(2,3) = (far+near)/(near-far);
    ortho(3,3) = 1;

    mglMultMatrix(ortho.values);
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
    color = vec4(red,green,blue,0);
}
