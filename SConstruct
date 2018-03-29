import os
env = Environment(ENV = os.environ)

env.Append(LIBS=["png"])
env.Append(CXXFLAGS=["-std=c++11","-g","-Wall","-O3"])
env.Append(LINKFLAGS=[])

env_nogl = env.Clone()
env_nogl.Append(CXXFLAGS=["-DNO_OPENGL"])

if env['PLATFORM'] == 'darwin':
    env.Append(CXXFLAGS=["-Wno-deprecated-declarations"])
    env.Append(LINKFLAGS=["-framework", "OpenGL", "-framework", "GLUT"])
else:
    env.Append(LIBS=["glut","GL","X11"])


minigl=env.Object("minigl.o","minigl.cpp");
dump_png=env.Object("dump_png.o","dump_png.cpp");
main=env.Object("main.o","main.cpp");
parse=env.Object("parse.o","parse.cpp");

env.Program("minigl",[main,minigl,parse,dump_png])

main_nogl=env_nogl.Object("main-nogl.o","main.cpp");
parse_nogl=env_nogl.Object("parse-nogl.o","parse.cpp");

env_nogl.Program("minigl-nogl",[main_nogl,minigl,parse_nogl,dump_png])
