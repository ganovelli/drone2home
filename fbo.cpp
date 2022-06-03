#include "fbo.h"
#include <GL/glew.h>
#include <stdio.h>




void FBO::Check(int fboStatus)
{
    switch(fboStatus){
    case GL_FRAMEBUFFER_COMPLETE_EXT:break;
    case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:          printf("FBO Incomplete: Attachment");break;
    case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:  printf("FBO Incomplete: Missing Attachment");break;
    case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:          printf("FBO Incomplete: Dimensions");break;
    case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:             printf("FBO Incomplete: Formats");break;
    case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:         printf("FBO Incomplete: Draw Buffer");break;
    case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:         printf("FBO Incomplete: Read Buffer");break;
    case GL_FRAMEBUFFER_UNSUPPORTED_EXT:                    printf("FBO Unsupported");break;
    default:                                                printf("Undefined FBO error");break;
    }
}


void FBO::Create(int w_,int h_)
{
    if( (w==w_) && (h==h_))
        return;
    w = w_;
    h = h_;
    glActiveTexture(GL_TEXTURE0);
	
	glGenTextures(1, &this->id_tex);
	glBindTexture(GL_TEXTURE_2D,  this->id_tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_DECAL);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	
    glGenFramebuffersEXT(1, &this->id_fbo);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, this->id_fbo);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, this->id_tex, 0);

    glGenRenderbuffersEXT(1, &this->id_depth);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT,this->id_depth);
    glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT,  w, h);
    glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, this->id_depth);

    int status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    Check(status);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, this->id_fbo);
    status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    Check(status);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void FBO::Remove()
{
 //	glDeleteTextures(1, &COLORBUFFERcalib);
    glDeleteFramebuffersEXT(1, &this->id_fbo);
    glDeleteRenderbuffersEXT(1, &this->id_depth);
}
