all: $(patsubst %.cc,%.out,$(wildcard *.cc))

%.out: %.cc
	g++ -std=c++20 \
	-DQUILE_ENABLE_LOGGING \
	-O3 -Wall -Wextra -pedantic \
	-I../../ -pthread \
	$< \
	&& ./a.out 2>&1 | fold >$@ \
	&& rm a.out \
	&& ./trim.sh $@
