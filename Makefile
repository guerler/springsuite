all:
	mkdir -p build
	g++ -O3 -ffast-math -lm -o build/TMalign spring_package/tmalign/TMalign.cpp
	cc -O3 -o build/pulchra spring_package/pulchra/pulchra.c spring_package/pulchra/pulchra_data.c -lm