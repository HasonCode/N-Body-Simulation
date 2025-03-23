CPPFLAGS = $(shell sdl2-config --cflags) $(shell $(PKG_CONFIG) SDL2_image --cflags) $(EXTRA_CPPFLAGS)
LDLIBS = $(shell sdl2-config --libs) $(shell $(PKG_CONFIG) SDL2_image --libs) -lGLEW $(EXTRA_LDLIBS)
EXTRA_LDLIBS ?= -lGL
PKG_CONFIG ?= pkg-config
all: cube
cube: shader_utils.o
clean:
	rm -f *.o triangle cube
.PHONY: all clean