# Compiler to use
CXX = g++

# Compiler flags:
# -std=c++11: Use C++11 standard (required for modern features)
# -Wall: Enable all warnings (good practice)
# -O3: Optimize the code for speed (important for simulations)
CXXFLAGS = -std=c++11 -Wall -O3

# The name of the final executable
TARGET = simulator

# The default rule (what happens when you type 'make')
all: $(TARGET)

# Rule to link the program
$(TARGET): main.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) main.cpp

# Rule to clean up files (what happens when you type 'make clean')
clean:
	rm -f $(TARGET)