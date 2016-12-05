TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

#LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    vec3.cpp \
    galacticcluster.cpp \
    celestialbody.cpp \
    integrator.cpp

HEADERS += \
    vec3.h \
    galacticcluster.h \
    celestialbody.h \
    integrator.h

QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS_RELEASE -= -O2
