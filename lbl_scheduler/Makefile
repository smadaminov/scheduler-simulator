all: dev 

dev:
	g++ -I/var/services/homes/smadaminov/libs/boost_1_66_0/ main.cpp -std=c++0x -fopenmp -o main.x

release:
	g++ -I/var/services/homes/smadaminov/libs/boost_1_66_0/ main.cpp -std=c++0x -O3 -fopenmp -o main.x

activate:
	g++ -I/var/services/homes/smadaminov/libs/boost_1_66_0/ main.cpp -I./libs -std=c++0x -O3 -fopenmp -o main.x -DACTIVATE -DDEBUG_PRINT

sbu_sched:
	g++ -g -I/var/services/homes/smadaminov/libs/boost_1_66_0/ main.cpp -g -I./libs -std=c++0x -O3 -fopenmp -o main.x -DSBU_SCHEDULER -DDEBUG_PRINT

lb_sched:
	g++ -I/var/services/homes/smadaminov/libs/boost_1_66_0/ main.cpp -I./libs -std=c++0x -O3 -fopenmp -o main.x -DLB_SCHEDULER -DDEBUG_PRINT

debug:
	g++ -g -I/var/services/homes/smadaminov/libs/boost_1_66_0/ main.cpp -I./libs -std=c++0x -O3 -fopenmp -o main.x -DDEBUG_PRINT

draw:
	g++ -I/var/services/homes/smadaminov/libs/boost_1_66_0/ main.cpp -std=c++0x -O3 -fopenmp -o main.x -DDEBUG_PRINT -DDRAW

clean:
	rm -f main.x
