# directories
INCLUDEDIR = include
SRCDIR = src
OBJDIR = obj
EXECUTABLEDIR = executables

# compiler
CC = g++
CFLAGS = -c -g -Wall `root-config --cflags` -I$(INCLUDEDIR)

#linker
LINKER = g++
LDFLAGS = `root-config --libs` -lSpectrum

# the := expands the meaning of the expression in the variable assignment 
SOURCES := $(wildcard $(SRCDIR)/*.cc) # take all the .cc files in the src folder
OBJECTS := $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o) # in the SOURCES (variable content) what matches $(SRCDIR)/%.cc becomes $(OBJDIR)/%.o where the % is used to create an "entry to entry" correspondance
TARGETS_SOURCES := $(wildcard $(SRCDIR)/*.cpp)
TARGETS_OBJECTS := $(TARGETS_SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
TARGETS := $(TARGETS_SOURCES:$(SRCDIR)/%.cpp=$(EXECUTABLEDIR)/%)

all: $(TARGETS)


$(TARGETS): $(EXECUTABLEDIR)/%: $(OBJDIR)/%.o $(OBJECTS)
	@$(LINKER) $< $(OBJECTS) $(LDFLAGS) -o $@
	@echo "\tLinking "$@" complete"

$(TARGETS_OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@$(CC) $(CFLAGS) $< -o $@
	@echo "\tCompiled "$<" succesfully"

# $< is the current input file and $@ is the target of this action the @ at the beginning of the line is to not print out the line
$(OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.cc # create a object $(OBJDIR)/%.o from the file $(SRCDIR)/%.cc for the name matching $(OBJDIR)/%.o in the OBJECT variable
	@$(CC) $(CFLAGS) $< -o $@
	@echo "\tCompiled "$<" succesfully"

clean:
	rm -f $(TARGETS) $(wildcard $(OBJDIR)/*)
