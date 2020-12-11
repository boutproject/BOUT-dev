BOUT_TOP	?= ../..

TARGET=backtrace

SOURCEC		= $(TARGET).cxx

include $(BOUT_TOP)/make.config

$(TARGET).o: $(TARGET).cxx
	@echo "  Compiling " $<
	@$(CXX) $(BOUT_INCLUDE) $(BOUT_FLAGS) -c $< -g -o $@
