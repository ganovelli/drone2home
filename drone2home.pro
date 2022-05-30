
CONFIG += console c++11
CONFIG -= app_bundle
QT -= gui
QT += xml
QT += opengl

VCGLIBPATH =../../libs/vcglib
INCLUDEPATH += $$VCGLIBPATH
INCLUDEPATH += $$VCGLIBPATH/eigenlib
INCLUDEPATH += ../../libs/AntTweakBar/include

# /home/ganovell/Documents/devel/libs/AntTweakBar/lib
LIBS += /home/ganovell/Documents/devel/libs/AntTweakBar/lib/libAntTweakBar.a
SOURCES += $$VCGLIBPATH/wrap/ply/plylib.cpp \
        main.cpp
