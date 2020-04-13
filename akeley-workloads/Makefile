all: 3mm-build cholesky-build particles cnn

3mm-build:
	gcc -O3 -I polybench-c-4.2.1-beta/utilities/ -I polybench-c-4.2.1-beta/linear-algebra/kernels/atax -DPOLYBENCH_USE_RESTRICT=1 -DPOLYBENCH_USE_C99_PROTO=1 -DMEDIUM_DATASET=1 polybench-c-4.2.1-beta/utilities/polybench.c polybench-c-4.2.1-beta/linear-algebra/kernels/3mm/3mm.c -o 3mm -lm

cholesky-build:
	gcc -O3 -I polybench-c-4.2.1-beta/utilities/ -I polybench-c-4.2.1-beta/linear-algebra/kernels/atax -DPOLYBENCH_USE_RESTRICT=1 -DPOLYBENCH_USE_C99_PROTO=1 -DMEDIUM_DATASET=1 polybench-c-4.2.1-beta/utilities/polybench.c polybench-c-4.2.1-beta/linear-algebra/solvers/cholesky/cholesky.c -o cholesky -lm

cnn: cnn-src/cnn.h cnn-src/cnn.cpp cnn-src/main.cpp
	c++ -g -O3 -Wall -Wextra cnn-src/cnn.cpp cnn-src/main.cpp -o cnn -lm

particles: particles-src/main.cpp particles-src/bedrock_particles.h particles-src/bedrock_particles.cpp
	$(CXX) -g -O3 -std=c++14 -Wall -Wextra particles-src/main.cpp particles-src/bedrock_particles.cpp -o particles -lm
