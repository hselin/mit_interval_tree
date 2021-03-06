
CC	=	g++
# use the following cflags for old versions of g++ with sets.old.H
#OLD_CFLAGS	=	-g -Wall -pedantic -fguiding-decls
# use the following for new versions of g++/egcs with sets.H
CFLAGS	=	-g -pedantic
LINK	=	g++
#LFLAGS	=	${CFLAGS} -ldl -ltcl -lefence
#OBJECTS = interval_tree.o test.o RedBlackTree.o
LFLAGS	=	${CFLAGS}
OBJECTS = interval_tree.o test.o

.PHONY: test clean 

#test: test.exe
#	./test.exe tree_test.tcl

test: $(OBJECTS)
	${LINK} ${LFLAGS} $(OBJECTS) -o test

clean:
	rm -f test *.o *~ MANIFEST

%.o : %.cc
	$(CC) $(CFLAGS) -c $< -o $@


