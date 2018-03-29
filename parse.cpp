#include <cstdio>
#include <cstdlib>
#include <vector>
#include <functional>
#include <cstring>
#include "minigl.h"

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

extern int width,height;

extern std::vector<std::function<void()> > commands_gl,commands_mgl;

void Parse(const char* test_file)
{
    FILE* F = fopen(test_file,"r");
    if(!F)
    {
        printf("Failed to open file '%s'\n",test_file);
        exit(EXIT_FAILURE);
    }

    float f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15;
    int i0,i1;
    char buff[1000],c;

    while(fgets(buff, sizeof(buff), F))
    {
        if(buff[0]=='#') continue;
        if(buff[0]=='\n') continue;
        if(sscanf(buff,"size %d %d\n",&i0,&i1)==2)
        {
            width=i0;
            height=i1;
        }
        else if(sscanf(buff,"glBegin GL_TRIANGLES%c\n",&c)==1)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glBegin(GL_TRIANGLES);});
#endif
            commands_mgl.push_back([=](){mglBegin(MGL_TRIANGLES);});
        }
        else if(sscanf(buff,"glBegin GL_QUADS%c\n",&c)==1)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glBegin(GL_QUADS);});
#endif
            commands_mgl.push_back([=](){mglBegin(MGL_QUADS);});
        }
        else if(sscanf(buff,"glEnd%c\n",&c)==1)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glEnd();});
#endif
            commands_mgl.push_back([=](){mglEnd();});
        }
        else if(sscanf(buff,"glVertex2 %g %g\n",&f0,&f1)==2)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glVertex2f(f0,f1);});
#endif
            commands_mgl.push_back([=](){mglVertex2(f0,f1);});
        }
        else if(sscanf(buff,"glVertex3 %g %g %g\n",&f0,&f1,&f2)==3)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glVertex3f(f0,f1,f2);});
#endif
            commands_mgl.push_back([=](){mglVertex3(f0,f1,f2);});
        }
        else if(sscanf(buff,"glMatrixMode GL_MODELVIEW%c\n",&c)==1)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glMatrixMode(GL_MODELVIEW);});
#endif
            commands_mgl.push_back([=](){mglMatrixMode(MGL_MODELVIEW);});
        }
        else if(sscanf(buff,"glMatrixMode GL_PROJECTION%c\n",&c)==1)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glMatrixMode(GL_PROJECTION);});
#endif
            commands_mgl.push_back([=](){mglMatrixMode(MGL_PROJECTION);});
        }
        else if(sscanf(buff,"glPushMatrix%c\n",&c)==1)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glPushMatrix();});
#endif
            commands_mgl.push_back([=](){mglPushMatrix();});
        }
        else if(sscanf(buff,"glPopMatrix%c\n",&c)==1)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glPopMatrix();});
#endif
            commands_mgl.push_back([=](){mglPopMatrix();});
        }
        else if(sscanf(buff,"glLoadIdentity%c\n",&c)==1)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glLoadIdentity();});
#endif
            commands_mgl.push_back([=](){mglLoadIdentity();});
        }
        else if(sscanf(buff,"glLoadMatrix %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                &f0,&f1,&f2,&f3,&f4,&f5,&f6,&f7,&f8,&f9,&f10,&f11,&f12,&f13,&f14,&f15)==16)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){float f[16]={f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15};glLoadMatrixf(f);});
#endif
            commands_mgl.push_back([=](){float f[16]={f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15};mglLoadMatrix(f);});
        }
        else if(sscanf(buff,"glMultMatrix %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                &f0,&f1,&f2,&f3,&f4,&f5,&f6,&f7,&f8,&f9,&f10,&f11,&f12,&f13,&f14,&f15)==16)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){float f[16]={f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15};glMultMatrixf(f);});
#endif
            commands_mgl.push_back([=](){float f[16]={f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15};mglMultMatrix(f);});
        }
        else if(sscanf(buff,"glTranslate %g %g %g\n",&f0,&f1,&f2)==3)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glTranslatef(f0,f1,f2);});
#endif
            commands_mgl.push_back([=](){mglTranslate(f0,f1,f2);});
        }
        else if(sscanf(buff,"glRotate %g %g %g %g\n",&f0,&f1,&f2,&f3)==4)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glRotatef(f0,f1,f2,f3);});
#endif
            commands_mgl.push_back([=](){mglRotate(f0,f1,f2,f3);});
        }
        else if(sscanf(buff,"glScale %g %g %g\n",&f0,&f1,&f2)==3)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glScalef(f0,f1,f2);});
#endif
            commands_mgl.push_back([=](){mglScale(f0,f1,f2);});
        }
        else if(sscanf(buff,"glFrustum %g %g %g %g %g %g\n",&f0,&f1,&f2,&f3,&f4,&f5)==6)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glFrustum(f0,f1,f2,f3,f4,f5);});
#endif
            commands_mgl.push_back([=](){mglFrustum(f0,f1,f2,f3,f4,f5);});
        }
        else if(sscanf(buff,"glOrtho %g %g %g %g %g %g\n",&f0,&f1,&f2,&f3,&f4,&f5)==6)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glOrtho(f0,f1,f2,f3,f4,f5);});
#endif
            commands_mgl.push_back([=](){mglOrtho(f0,f1,f2,f3,f4,f5);});
        }
        else if(sscanf(buff,"glColor %g %g %g\n",&f0,&f1,&f2)==3)
        {
#ifndef NO_OPENGL
            commands_gl.push_back([=](){glColor3f(f0,f1,f2);});
#endif
            commands_mgl.push_back([=](){mglColor(f0,f1,f2);});
        }
        else
        {
            int len=strlen(buff);
            if(buff[len-1]=='\n') buff[len-1]=0;
            printf("Unrecognized command: '%s'\n",buff);
            exit(EXIT_FAILURE);
        }
    }
    fclose(F);
}
