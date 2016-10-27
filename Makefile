
# File Locations
SOURCES_ROOT:=./
BINDIR:=./

# File Information
SOURCES_EXT=c
HEADER_EXT=h
EXEC=sm_parallel

# JCL Control filename
JCL=./sm_job.slm

# Include Directories
INCLUDE_DIRS=\
.\

# Compiler Information
CC:=gcc
CSTD:=-std=c99
CWARN:=-Wall -Wextra -Werror
DEBUG_FLAGS:=-g -fsanitize=thread -fsanitize=undefined
RELEASE_FLAGS:=-O3

# Compiler Flags Generation
CFLAGS:=$(CSTD) $(CWARN)

# Linker Flags
LDFLAGS:=

# Executable Setup
CLEAN=rm -f

# Libraries
LIBS:= \

PTHREAD:=-pthread

# File Lists
SOURCES:=$(wildcard $(SOURCES_ROOT)/*.$(SOURCES_EXT))
HEADERS:=$(wildcard $(SOURCES_ROOT)/*.$(HEADER_EXT))

# Objects List
OBJS:=$(SOURCES:%.c=%.o)

# Make Targets
.PHONY: default run release batch clean

default: CFLAGS += $(DEBUG_FLAGS)
default: LDFLAGS += $(DEBUG_FLAGS)
default: $(OBJS)
default: mkbin
	@echo "Building with Debugging"
	$(CC) $(OBJS) $(LDFLAGS) -o $(BINDIR)$(EXEC) $(LIBS:%=-l%) $(PTHREAD)

run: default
	@echo "Running"
	@if [ -f $(BINDIR)/$(EXEC) ]; then \
		$(BINDIR)/$(EXEC); \
	else \
		echo "Bad build, no run target found."; \
	fi

mkbin:
	@if [ ! -d "$(BINDIR)" ]; then \
		mkdir $(BINDIR); \
	fi

release: CFLAGS += $(RELEASE_FLAGS)
release: LDFLAGS += $(RELEASE_FLAGS)
release: $(OBJS)
release: mkbin
	@echo "Building Release"
	$(CC) $(OBJS) $(LDFLAGS) -o $(BINDIR)$(EXEC) $(LIBS:%=-l%) $(PTHREAD)

batch: release
	@echo "Submitting Batch Job to Balena"
	sbatch $(BATCH)

clean:
	@echo "Cleaning build artefacts"
	$(RM) $(OBJS:%=$(SOURCES_ROOT)/%)
	$(RM) $(BINDIR)/$(EXEC)

# Compile Sources
$(SOURCES_ROOT)/%.o: $(SOURCES_ROOT)%.$(SOURCES_EXT)
	$(CC) $(CFLAGS) $(INCLUDE_DIRS:%=-I%) $(PTHREAD) -c -o $@ $<

