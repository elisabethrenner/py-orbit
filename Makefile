DIRS  = 
DIRS += ./ext
DIRS += ./src

#========rules=========================

all: compile

ifeq ($(COMP),)
    COMP=intel
endif

$(info COMP is $(COMP))

compile:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then echo "Compiling $$dir"; $(MAKE) compile; fi;); \
	done	

clean:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) clean; fi;); \
	done
	rm -rf ./bin/pyORBIT
	rm -rf ./doc/html

docs:
	doxygen
