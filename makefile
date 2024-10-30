# makefile

INCLUDES = -I./
DEFINES = 

CC := gcc
LD := gcc
LINKFLAGS := 
CFLAGS := -Wall -Wextra -Iinclude -fdiagnostics-color=always -std=c17 $(INCLUDES) $(DEFINES) -lm -MMD -MP

BUILD_DIR := build

DEBUG_LIBS := -lm
DEBUG_CFLAGS := $(CFLAGS) -g
DEBUG_LINKFLAGS := $(LINKFLAGS) -g
DEBUG_DIR := $(BUILD_DIR)/debug

RELEASE_LIBS := -lm
RELEASE_CFLAGS := $(CFLAGS) -O3
RELEASE_LINKFLAGS := $(LINKFLAGS)
RELEASE_DIR := $(BUILD_DIR)/release

# Make does not offer a recursive wildcard function, so here's one:
# DO NOT ADD SPACES TO ARGS WHILE CALLING
# leave 2nd arg empty to search from this directory
rwildcard=$(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))

EXCLUDE := test/%.c

DEBUG_SRC_C := $(filter-out $(EXCLUDE),$(call rwildcard,,*.c))
DEBUG_OBJ_C := $(patsubst %.c, $(DEBUG_DIR)/objs/%.o, $(DEBUG_SRC_C))
DEBUG_DEPENDS := $(patsubst %.c,$(DEBUG_DIR)/objs/%.d,$(DEBUG_SRC_C))

RELEASE_SRC_C := $(filter-out $(EXCLUDE),$(call rwildcard,,*.c))
RELEASE_OBJ_C := $(patsubst %.c, $(RELEASE_DIR)/objs/%.o, $(RELEASE_SRC_C))
RELEASE_DEPENDS := $(patsubst %.c,$(RELEASE_DIR)/objs/%.d,$(RELEASE_SRC_C))

.PHONY: all debug release test clean always

all: build_debug

debug: build_debug
release: build_release

build_debug: $(DEBUG_DIR)/prog
build_release: $(RELEASE_DIR)/prog

$(DEBUG_DIR)/prog: $(DEBUG_OBJ_C)
	@echo "Linking Project(Debug)"
	@$(LD) $(DEBUG_LINKFLAGS) -o $@ $^ $(DEBUG_LIBS)
	@echo "Created " $@
	@cp $@ $(BUILD_DIR)/prog

-include $(DEBUG_DEPENDS)

$(DEBUG_DIR)/objs/%.o: %.c makefile
	@mkdir -p $(@D)
	@echo "[Debug] Compiling: " $<
	@$(CC) $(DEBUG_CFLAGS) -c -o $@ $<


$(RELEASE_DIR)/prog: $(RELEASE_OBJ_C)
	@echo "Linking Project(Release)"
	@$(LD) $(RELEASE_LINKFLAGS) -o $@ $^ $(RELEASE_LIBS)
	@echo "Created " $@

-include $(RELEASE_DEPENDS)

$(RELEASE_DIR)/objs/%.o: %.c makefile
	@mkdir -p $(@D)
	@echo "[Release] Compiling: " $<
	@$(CC) $(RELEASE_CFLAGS) -c -o $@ $<

clean:
	@rm -r $(BUILD_DIR)


