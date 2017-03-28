LIB_SRC = ../lib/bbpr.c
LIB_H   = ../lib/

trapezoid: ./trapezoid/main.c ./trapezoid/motionPlanning.c
	gcc $^ $(LIB_SRC) -I$(LIB_H) -I./trapzeoid/ -o test 

s-curve: ./s-curve/main.c ./s-curve/motionPlanning.c
	gcc $^ $(LIB_SRC) -I$(LIB_H) -I./s-curve/ -o test 

5th-order: ./5th-order/main.c ./5th-order/motionPlanning.c
	gcc $^ $(LIB_SRC) -I$(LIB_H) -I./5th-order/ -o test 

plot:
	./test > test.out && gnuplot plot.gnu && open plot.png