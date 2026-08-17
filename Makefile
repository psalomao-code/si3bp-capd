CXX ?= g++
CAPD_CONFIG ?= capd-config
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic

all: si3bp_capd

si3bp_capd: si3bp_capd.cpp
	$(CXX) $(CXXFLAGS) $< $$($(CAPD_CONFIG) --cflags --libs) -o $@

run: si3bp_capd
	./si3bp_capd | tee capd-proof.log

clean:
	rm -f si3bp_capd capd-proof.log
