# Makefile for Geode
#
# Copyright (C) 2003, 2006 Jesse H. Willett
# email: jhw@lmi.net, jhwjhw@gmail.com
# web:   http://users.lmi.net/~jhw

SHELL=/bin/sh
.SUFFIXES:

INCLUDE := -I/usr/X11R6/include/
LIBDIR  := -L/usr/X11R6/lib -L/usr/lib
LIBS += -lX11 
LIBS += -lXi 
LIBS += -lXmu 
LIBS += -lm 
LIBS += -lGL 
LIBS += -lGLU 
LIBS += -lglut 
LIBS += -lstdc++ 
LIBS += -lsupc++
ARGS += -g 
ARGS += -pedantic 
ARGS += -Wall
#ARGS += -pg

SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.hpp)
OBJECTS = $(SOURCES:%.cpp=build/%.o)

TARGET := bin/geodefunc

.PHONY: all
all: $(TARGET)

.PHONY: clean
clean: cl
	@echo cleaning build targets
	rm -rf $(TARGET)
	rm -rf build

.PHONY: cl
cl:
	@echo massaging out editing artifacts
	rm -rf core* *~ 
	chmod -xxx *.hpp *.cpp Makefile
	chmod +xxx bin/*

$(TARGET): $(OBJECTS)
	@echo linking "$(dir $<)*.o" to $@
	@mkdir -p $(dir $@)
	@$(CC) $(ARGS) $(LIBDIR) $^ -o $@ $(LIBS)

$(OBJECTS): build/%.o: %.cpp $(HEADERS)
	@echo $<" --$(CC)--> "$@
	@mkdir -p $(dir $@)
	@$(CC) $(ARGS) $(INCLUDE) -c $< -o $@
