# Makefile for Geode
#
# Copyright (C) 2003, 2006 Jesse H. Willett
# email: jhw@lmi.net, jhwjhw@gmail.com
# web:   http://users.lmi.net/~jhw

SHELL=/usr/bin/sh
.SUFFIXES:

INCLUDE = -I/usr/X11R6/include/
LIBDIR  = -L/usr/X11R6/lib -L/usr/lib
LIBS = -lX11 -lXi -lXmu -lm -lGL -lGLU -lglut -lstdc++ -lsupc++
ARGS = -g -pedantic -Wall #-pg

SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.hpp)
OBJECTS = $(SOURCES:%.cpp=%.o)

TARGET = bin/geodefunc

all: $(TARGET)

clean:
	rm -f *.o core* *~ $(TARGET) gmon.out

$(TARGET): $(OBJECTS)
	$(CC) $(ARGS) $(LIBDIR) $^ -o $@ $(LIBS)

$(OBJECTS): %.o: %.cpp $(HEADERS)
	$(CC) $(ARGS) $(INCLUDE) -c $< -o $@
