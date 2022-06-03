#ifndef FBO_UTIL
#define FBO_UTIL
#include <GL/glew.h>
#include <iostream>
struct FBO{
    int w,h;
    GLuint id_fbo,id_tex,id_depth;
    void Check(int fboStatus);
    void Create(int,int);
    void Remove();
};

#endif
