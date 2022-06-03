#ifndef _SHADER_BASIC_
#define _SHADER_BASIC_

#include <GL/glew.h>
#include <stdio.h>
#include <string>
#include <map>

struct Shader{
        GLuint   vs, gs, fs, pr;

        std::map<std::string,int> uni;
        void bind(std::string name){
            uni[name] = glGetUniformLocation(pr, name.c_str());
        }

        int operator[](std::string name){
            return uni[name];
        }

        static  char *textFileRead(const char *fn) {

		FILE *fp;
		char *content = NULL;

		int count = 0;

		if (fn != NULL) {
			fp = fopen(fn, "rt");

			if (fp != NULL) {

				fseek(fp, 0, SEEK_END);
				count = ftell(fp);
				rewind(fp);

				if (count > 0) {
					content = (char *)malloc(sizeof(char)* (count + 1));
					count = fread(content, sizeof(char), count, fp);
					content[count] = '\0';
				}
				fclose(fp);
			}
		}
		return content;
	}

        static int  SetFromString( const GLchar *stringV,const GLchar *stringG,  const GLchar *stringF,  GLuint & vs,GLuint &gs,GLuint &fs, GLuint &pr){
                fs = glCreateShader(GL_FRAGMENT_SHADER);
                vs = glCreateShader(GL_VERTEX_SHADER);
                gs = glCreateShader(GL_GEOMETRY_SHADER);

		if(stringV!=NULL){
			glShaderSource(vs, 1, &stringV,NULL);
			glCompileShader(vs);
			int errV;
			glGetShaderiv(vs,GL_COMPILE_STATUS,&errV);
			if(errV!=GL_TRUE) return -2;
		}else return -1;

                if(stringG!=NULL){
                        glShaderSource(gs, 1, &stringG,NULL);
                        glCompileShader(gs);
                        int errG;
                        glGetShaderiv(gs,GL_COMPILE_STATUS,&errG);
                        if(errG!=GL_TRUE) return -5;
                };

		if(stringF!=NULL){
			glShaderSource(fs, 1, &stringF,NULL);
			glCompileShader(fs);
			int errF;
			glGetShaderiv(fs,GL_COMPILE_STATUS,&errF);
			if(errF!=GL_TRUE) return -4;
		}else return -3;
		pr = glCreateProgram();
		glAttachShader(pr,vs);
		glAttachShader(pr,fs);
                if(stringG!=NULL)
                    glAttachShader(pr,gs);

		glLinkProgram(pr);
		return 0;
	}

        int  SetFromFile( const GLchar *nameV, const char *nameG, const char *nameF){
		fs= glCreateShader(GL_FRAGMENT_SHADER);
		vs= glCreateShader(GL_VERTEX_SHADER);

                char * code;

		if(nameV!=NULL){
			code = textFileRead(nameV);
			glShaderSource(vs, 1, &code,NULL);
			glCompileShader(vs);
			int n;
			char b[65536];
			glGetShaderInfoLog(vs, 65536, &n, b);
			b[n] = '\0';

			int errV;
			glGetShaderiv(vs,GL_COMPILE_STATUS,&errV);
			if(errV!=GL_TRUE)
				return -2;
		}else return -1;

                if(nameG!=NULL){
                         gs = glCreateShader(GL_GEOMETRY_SHADER);
                        code = textFileRead(nameG);
                        glShaderSource(gs, 1, &code,NULL);
                        glCompileShader(gs);
                        int n;
                        char b[65536];
                        glGetShaderInfoLog(gs, 65536, &n, b);
                        b[n] = '\0';

                        int errG;
                        glGetShaderiv(vs,GL_COMPILE_STATUS,&errG);
                        if(errG!=GL_TRUE)
                                return -5;
                };

		if(nameF!=NULL){
			code = textFileRead(nameF);
			glShaderSource(fs, 1, &code,NULL);
			glCompileShader(fs);
			int n;
			char b[65536];
			glGetShaderInfoLog(fs, 65536, &n, b);
			b[n] = '\0';

			int errF;
			glGetShaderiv(fs,GL_COMPILE_STATUS,&errF);
			if(errF!=GL_TRUE) 
				return -4;
		}else return -3;
		pr = glCreateProgram();
                if(nameV!=NULL) glAttachShader(pr,vs);
                if(nameG!=NULL) glAttachShader(pr,gs);
                if(nameF!=NULL) glAttachShader(pr,fs);
		glLinkProgram(pr);
                glValidateProgram(pr);
                GLuint err = 0;
                err =  glGetError();
                char res[65536];
                int l;
                glGetProgramInfoLog(pr,65536,&l,res);
                res[l] =  '\0';
		return 0;
	}

        void Validate( ){
		int res;
                const GLuint & s = this->pr;
		glValidateProgram(s);
		glGetProgramiv(s,GL_VALIDATE_STATUS,&res);
		printf("validation of program %d:%d \n",s,res);

		glGetProgramiv(s,GL_LINK_STATUS,&res);
		printf("linking of program %d:%d \n",s,res);

		glGetProgramiv(s,GL_ACTIVE_ATTRIBUTES,&res);
		printf("active attribute of program %d:%d \n",s,res);

		glGetProgramiv(s,GL_ACTIVE_UNIFORMS,&res);
		printf("active uniform  of program %d:%d \n",s,res);

		glGetProgramiv(s,GL_ACTIVE_UNIFORM_MAX_LENGTH,&res);
		printf("active uniform Max Length of program %d:%d \n",s,res);
}
};

#endif
