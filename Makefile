CC = gcc
SOURCEDIR = src/
BUILDDIR = build

EXECUTABLE = scalaf
SOURCES = $(wildcard src/*.c)
OBJECTS = $(patsubst $(SOURCEDIR)/%.c,$(BUILDDIR)/%.o, $(SOURCES))

LDFLAGS = -lm

all: dir $(BUILDDIR)/$(EXECUTABLE)

dir: 
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/$(EXECUTABLE): $(OBJECTS) 
	$(CC) $^ -o $@ $(LDFLAGS)

clean: 
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/$(EXECUTABLE)