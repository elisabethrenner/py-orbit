DIRS  = 
#DIRS += ./ecloud
#DIRS += ./intsystem
#DIRS += ./laserstripping
#DIRS += ./H0starkeffect
DIRS += ./PTC

$(info Got $$LIBS as follows [${LIBS}])

#========rules=========================

all: compile

compile:  
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile.$(COMP) ]; then echo "Compiling $$dir"; $(MAKE) -f Makefile.$(COMP) compile; fi;); \
	done
	@echo "Compilation of py-orbit/ext sub-directories completed"
	@echo "Compilation of py-orbit/src completed"
	@echo "****************************************"

clean:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) clean; fi;); \
	done
