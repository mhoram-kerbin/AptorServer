TARGET := AptorServer

MODULE_DIR = Optimizer
MODULES := $(wildcard $(MODULE_DIR)/*.cpp)
MODULES_OBJS = $(MODULES:%.cpp=%.o)
MODULES_HEADER = $(MODULES:%.cpp=%.h)

INTERFACES_DIR = Interface
INTERFACES := $(wildcard $(INTERFACES_DIR)/*.cpp)
INTERFACES_OBJS = $(INTERFACES:%.cpp=%.o)
INTERFACES_HEADER = $(INTERFACES:%.cpp=%.h)

EXTRALIBS = -lModernPsoptInterface -lws2_32

PARTS := Optimizer.o Control.o

OBJECT_DEPS = $(PARTS) $(MODULES_OBJS) $(INTERFACES_OBJS)

EXTRAHEADER = Exception/ServerException.h

CXXFLAGSEXTRA = -std=c++0x -pedantic -iquote.

include $(PSOPT_DIR)/Makefile_include.mk

re: psoptclean $(TARGET)

Optimizer.o: Optimizer.cpp Optimizer.h $(EXTRAHEADER)
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGSEXTRA) $(INCLUDES) $< -o $@

Control.o: Control.cpp Control.h $(EXTRAHEADER)
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGSEXTRA) $(INCLUDES) $< -o $@

$(MODULE_DIR)/%.o: $(MODULE_DIR)/%.cpp $(MODULE_DIR)/%.h Optimizer.h $(EXTRAHEADER)
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGSEXTRA) $(INCLUDES) $< -o $@

$(INTERFACES_DIR)/%.o: $(INTERFACES_DIR)/%.cpp $(INTERFACES_DIR)/%.h $(EXTRAHEADER)
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGSEXTRA) $(INCLUDES) $< -o $@

clean2:
	rm -f *.o $(TARGET) $(TARGET).exe $(APP_OBJS)

projectclean:
	rm -f $(MODULES_OBJS) $(INTERFACES_OBJS)

projectpsoptclean:
	rm -f *.pdf *.png
