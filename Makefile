all:
	mkdir -p build
	g++ -O3 -ffast-math -lm -o build/TMalign dependencies/tmalign/TMalign.cpp
	cc -O3 -o build/pulchra dependencies/pulchra/pulchra.c dependencies/pulchra/pulchra_data.c -lm
