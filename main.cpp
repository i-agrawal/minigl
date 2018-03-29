/**
 * main.cpp
 * -------------------------------
 * This is simple testbed for your miniGL implementation.
 *
 * Two versions of the program are produced.  One compiles with OpenGL (minigl)
 * and the other does not (minigl-nogl).  The latter is useful when running
 * over SSH.  It will also be used for grading.
 *
 * Usage: minigl -i <input-file> [ -s <solution-file> ] [ -o <stats-file> ] [ -c ] [ -p ] [ -g ]
 *     <input-file>      File with commands to run
 *     <solution-file>   File with solution to compare with
 *     <stats-file>      Dump statistics to this file rather than stdout
 *     -p                Dump results to png files
 *     -g                Use OpenGL
 *     -c                Compare only; if using OpenGL, don't leave window open
 *
 * Only the -i is manditory.  You must specify a test to run.  For example:
 *
 * minigl -i 00.txt
 *
 * Some potentially useful ways to run the program are described below.
 *
 * Compare visually with OpenGL solution
 *
 * minigl -gi 00.txt
 *
 * Compare visually with OpenGL solution; dump your solution, OpenGL solution,
 * and the difference between the two to png files.  This is useful since it
 * lets you examine the error visually.
 *
 * minigl -gpi 00.txt
 *
 * Same, but do not keep the OpenGL window open (when you only need the png)
 *
 * minigl -cgpi 00.txt
 *
 * If you want to test over SSH, you can use mingl-nogl, which does not have
 * any OpenGL dependence.  The solution can be specified as a file.  For
 * example, this will compare your solution to 00.png, outputing your
 * solution and the error to png.  When the solution is specified as a png,
 * it is not dumped back out as png.
 *
 * minigl-nogl -i 00.txt -s 00.png -p
 *
 * The -o options is primarily intended for grading.  It allows the statistics
 * to be accessed without being affected by debugging output.  The -n option
 * runs the commands multiple times to get a better timing estimate.
 *
 * Grading will be done like this:
 *
 * minigl-nogl -i 00.txt -s 00.png -o result.txt -n 3
 *
 * Your code is timed each time it is run.  The output is reported in
 * milliseconds.  This is how your solution's runtime will be determined.
 */
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include "minigl.h"
#include <functional>
#include <cassert>
#include <climits>
#include <unistd.h>

#ifndef NO_OPENGL
#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif
#endif

#ifdef __MACH__
#include <sys/time.h>
//clock_gettime is not implemented on OSX
#define CLOCK_MONOTONIC 0
int clock_gettime(int /*clk_id*/, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#endif

// Change these to test MiniGL at different resolutions
int width=320,height=240;
MGLpixel* pixel_data;
MGLpixel* pixel_data_sol;
bool show_sol = false;

std::vector<std::function<void()> > commands_gl,commands_mgl;

void Parse(const char* test_file);
void Dump_png(MGLpixel* data,int width,int height,const char* filename);
void Read_png(MGLpixel*& data,int& width,int& height,const char* filename);

#ifndef NO_OPENGL
void Display_Side_By_Side()
{
    glViewport(0,0,width,height);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    // Generate solution with OpenGL
    if(!show_sol)
    {
        for(size_t i=0;i<commands_gl.size();i++)
            commands_gl[i]();
    }
    
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if(show_sol)
    {
        glRasterPos2f(0,0);
        glDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_INT_8_8_8_8,(GLvoid*)pixel_data_sol);
    }

    glRasterPos2f(1,0);
    glDrawPixels(width,height,GL_RGBA,GL_UNSIGNED_INT_8_8_8_8,(GLvoid*)pixel_data);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glFlush();
}

void Capture_OpenGL_Solution()
{
    glReadPixels(0,0,width,height,GL_RGBA,GL_UNSIGNED_INT_8_8_8_8,(GLvoid*)pixel_data_sol);
}
#endif

void Compare(FILE* stats_file, bool use_png)
{
    std::string file_gl="gl.png",file_mgl="mgl.png",file_diff="diff.png";
    MGLpixel pixel_data_diff[width*height];
    for(int i=0;i<width*height;i++) pixel_data_diff[i]=Make_Pixel(0,0,0);
    int total_diff = 0;
    for(int i=0;i<width*height;i++)
    {
        int A=pixel_data_sol[i];
        int B=pixel_data[i];
        int rA,gA,bA,rB,gB,bB;
        From_Pixel(A,rA,gA,bA);
        From_Pixel(B,rB,gB,bB);
        int r=abs(rA-rB);
        int g=abs(gA-gB);
        int b=abs(bA-bB);
        int diff=std::max(std::max(r,g),b);
        total_diff += diff;
        MGLpixel diff_color=Make_Pixel(diff, diff, diff);
        pixel_data_diff[i]=diff_color;
    }
    if(use_png) Dump_png(pixel_data_diff,width,height,"diff.png");
    fprintf(stats_file, "diff: %.2f\n",total_diff/(2.55*width*height));
}

void Usage(const char* prog_name)
{
    std::cerr<<"Usage: "<<prog_name<<" -i <input-file> [ -s <solution-file> ] [ -o <stats-file> ] [ -c ] [ -p ] [ -g ]"<<std::endl;
    std::cerr<<"    <input-file>      File with commands to run"<<std::endl;
    std::cerr<<"    <solution-file>   File with solution to compare with"<<std::endl;
    std::cerr<<"    <stats-file>      Dump statistics to this file rather than stdout"<<std::endl;
    std::cerr<<"    -p                Dump results to png files"<<std::endl;
    std::cerr<<"    -g                Use OpenGL"<<std::endl;
    std::cerr<<"    -c                Compare only; if using OpenGL, don't leave window open"<<std::endl;
    std::cerr<<"    -n                Run solution multiple times to get more reliable timing results"<<std::endl;
    exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
    const char* solution_file = 0;
    const char* input_file = 0;
    const char* statistics_file = 0;
    bool use_png = false;
    bool use_opengl = false;
    bool compare_only = false;
    int repeat=1;
    
    // Parse commandline options
    while(1)
    {
        int opt = getopt(argc, argv, "s:i:pgo:cn:");
        if(opt==-1) break;
        switch(opt)
        {
            case 's': solution_file = optarg; break;
            case 'i': input_file = optarg; break;
            case 'p': use_png = true; break;
#ifndef NO_OPENGL
            case 'g': use_opengl = true; break;
#endif
            case 'o': statistics_file = optarg; break;
            case 'c': compare_only = true; break;
            case 'n': repeat = atoi(optarg); break;
        }
    }

    if(compare_only && solution_file) use_opengl = false;

    // Sanity checks
    if(!input_file)
    {
        std::cerr<<"Test file required.  Use -i."<<std::endl;
        Usage(argv[0]);
    }

    // Parse the OpenGl commands
    Parse(input_file);

    FILE* stats_file = stdout;
    if(statistics_file) stats_file = fopen(statistics_file, "w");

    // Run and time the student solution
    pixel_data = new MGLpixel[width*height];

    long long time_diff = LLONG_MAX;
    for(int j=0;j<repeat;j++)
    {
        for(int i=0;i<width*height;i++) pixel_data[i]=Make_Pixel(0,0,0);

        timespec t0,t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        for(size_t i=0;i<commands_mgl.size();i++)
            commands_mgl[i]();
        mglReadPixels(width,height,pixel_data);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        long long a0 = (long long)t0.tv_sec*1000000000+t0.tv_nsec;
        long long a1 = (long long)t1.tv_sec*1000000000+t1.tv_nsec;
        time_diff = std::min(time_diff, a1-a0);
    }
    fprintf(stats_file, "time: %.2f\n",time_diff*1e-6);

    // Solution
    if(solution_file || use_opengl)
    {
        pixel_data_sol = new MGLpixel[width*height];
        for(int i=0;i<width*height;i++) pixel_data_sol[i]=Make_Pixel(0,0,0);
    }

    // Read solution from file
    if(solution_file)
    {
        int width_sol = 0;
        int height_sol = 0;
        Read_png(pixel_data_sol, width_sol, height_sol, solution_file);
        if(width!=width_sol || height!=height_sol)
        {
            std::cerr<<"Solution dimensions ("<<width_sol<<","<<height_sol<<") do not match problem size ("<<width<<","<<height<<")"<<std::endl;
            exit(EXIT_FAILURE);
        }
        show_sol = true;
    }

#ifndef NO_OPENGL
    if(use_opengl)
    {
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);
        glutInitWindowSize(2*width, height);
        glutCreateWindow("Test OpenGL Program");
        glClearColor(0,0,0,1);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glShadeModel(GL_SMOOTH);
        glPushAttrib(GL_COLOR_BUFFER_BIT);
        glDrawBuffer(GL_FRONT);

        Display_Side_By_Side();
        if(!show_sol) Capture_OpenGL_Solution();
    }
#endif

    if(use_png)
    {
        Dump_png(pixel_data,width,height,"mgl.png");
        if(use_opengl && !solution_file)
            Dump_png(pixel_data_sol,width,height,"gl.png");
    }
    
    // Compare solutions
    if(pixel_data_sol) Compare(stats_file, use_png);
    
#ifndef NO_OPENGL
    // Leave results up for display.
    if(use_opengl && !compare_only)
    {
        glutDisplayFunc(Display_Side_By_Side);
        glutMainLoop();
    }
#endif

    if(stats_file != stdout) fclose(stats_file);
    delete [] pixel_data;
    delete [] pixel_data_sol;
    return 0;
}
