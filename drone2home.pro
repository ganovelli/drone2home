
CONFIG += console c++11
CONFIG -= app_bundle
QT -= gui
QT += xml
QT += opengl

VCGLIBPATH =../../libs/vcglib
INCLUDEPATH += $$VCGLIBPATH
INCLUDEPATH += $$VCGLIBPATH/eigenlib

SOURCES += $$VCGLIBPATH/wrap/ply/plylib.cpp \
        main.cpp
