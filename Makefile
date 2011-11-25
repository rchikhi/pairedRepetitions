OBJS = \
       my_sarray/lcp.o \
       my_sarray/sarray.o \
       my_sarray/scode.o 
LIBS = -L ./judy-1.0.5/lib -L ./argtable2-11/lib
LFLAGS = -lJudy -largtable2
MY_FILES_V8 = 33bits.c byte2judy.c rb-1bit.c
MY_FILES = 33bits.c short2judy.c rb-1bit.c


.PHONY:run clean

default: pairs-v7
	./pairs-v7

v8: pairs-v8
	./pairs-v8

v8mp: pairs-v8mp

array: pairs-v7-array
	./pairs-v7

multi: pairs-v7-multi
	./pairs-v7

2bits: pairs-v7-2bits
	./pairs-v7

simple: pairs-v7-simple
	./pairs-v7

memcheck: $(OBJS)
	rm -f callgrind*
	valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./pairs-v7

memprofile: $(OBJS)
	rm -f callgrind*
	valgrind --tool=massif --massif-out-file=massif.out.here  ./pairs-v7


callprofile: $(OBJS)
	rm -f callgrind*
	valgrind --callgrind-out-file=callgrind.out.here --tool=callgrind --dump-instr=yes ./pairs-v7-hd

pairs-v7: $(OBJSj) pairs-v7.c p-judy.c
	rm -f *.o pairs-v7
	gcc  -O3 -static -I ./judy-1.0.5/include $(LIBS) -o pairs-v7 $(OBJS) pairs-v7.c p-judy.c $(MY_FILES) $(LFLAGS)

pairs-v8: $(OBJSj) pairs-v8.c p-judy-seq.c
	rm -f *.o pairs-v8
	gcc43  -O3 -static -I ./judy-1.0.5/include $(LIBS) -o pairs-v8 $(OBJS) pairs-v8.c p-judy-seq.c $(MY_FILES) $(LFLAGS) -Wall

pairs-v8mp: $(OBJSj) pairs-v8mp.c p-judy.c
	rm -f *.o pairs-v8mp
	gcc43 -fopenmp  -O3 -static -I ./judy-1.0.5/include $(LIBS) -o pairs-v8mp $(OBJS) pairs-v8mp.c p-judy.c $(MY_FILES) $(LFLAGS)


pairs-v7-array: $(OBJSj) pairs-v7.c p-array.c
	rm -f *.o pairs-v7
	gcc  -O3 -static -I ./judy-1.0.5/include $(LIBS) -o pairs-v7 $(OBJS) pairs-v7.c p-array.c $(MY_FILES) $(LFLAGS)


pairs-v7-multi: $(OBJSj) pairs-v7.c p-judy-multi.c
	rm -f *.o pairs-v7
	gcc  -O3 -static -I ./judy-1.0.5/include -L ./judy-1.0.5/lib -o pairs-v7 $(OBJS) pairs-v7.c p-judy-multi.c  $(MY_FILES) -lJudy

pairs-v7-2bits: $(OBJSj) pairs-v7.c p-judy-2bits-array.c
	rm -f *.o pairs-v7
	gcc  -O3 -static -I ./judy-1.0.5/include $(LIBS) -o pairs-v7 $(OBJS) pairs-v7.c p-judy-2bits-array.c $(MY_FILES) $(LFLAGS)


pairs-v7-simple: $(OBJSj) pairs-v7.c p-judy-simple.c
	rm -f *.o pairs-v7
	gcc  -O3 -static -I ./judy-1.0.5/include -L ./judy-1.0.5/lib -o pairs-v7 $(OBJS) pairs-v7.c p-judy-simple.c $(MY_FILES) -lJudy




clean:
	rm -f *.o pairs-v7 pairs-v7-hd
