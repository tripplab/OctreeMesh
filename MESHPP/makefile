CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic

COMMON_SRC = src/post_msh_reader.cc src/post_msh_writer.cc src/operations.cc
ROUNDTRIP_SRC = src/roundtrip_main.cc
APPLY_SRC = src/apply_main.cc

all: meshpp_roundtrip meshpp_apply

meshpp_roundtrip: $(COMMON_SRC:.cc=.o) $(ROUNDTRIP_SRC:.cc=.o)
	$(CXX) $(CXXFLAGS) -o $@ $^

meshpp_apply: $(COMMON_SRC:.cc=.o) $(APPLY_SRC:.cc=.o)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f src/*.o meshpp_roundtrip meshpp_apply

SAN_FLAGS = -fsanitize=address,undefined -fno-omit-frame-pointer

sanitize:
	$(CXX) $(CXXFLAGS) $(SAN_FLAGS) -o meshpp_roundtrip_asan src/post_msh_reader.cc src/post_msh_writer.cc src/operations.cc src/roundtrip_main.cc
	$(CXX) $(CXXFLAGS) $(SAN_FLAGS) -o meshpp_apply_asan src/post_msh_reader.cc src/post_msh_writer.cc src/operations.cc src/apply_main.cc

clean-sanitize:
	rm -f meshpp_roundtrip_asan meshpp_apply_asan
