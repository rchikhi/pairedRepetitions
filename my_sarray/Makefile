OBJS = \
       lcp.o \
       sarray.o \
       scode.o \
       test.o

.PHONY:runtest clean
runtest:test
	./test

test:$(OBJS)
	gcc -o test $(OBJS)

clean:
	rm *.o test
