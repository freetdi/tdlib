CC=g++
SOFLAGS=-shared -Wl,-soname,libtd.so


all: interface cythonize pylib partial_clean test

interface: python_tdlib.cpp python_tdlib.hpp
	$(CC) $(LDFLAGS) $(SOFLAGS) -I./../lib -fPIC -o libtd.so python_tdlib.cpp -c

cythonize: tdlib.pyx tdlib.pxd
	cython --cplus tdlib.pyx
	$(CC) -pthread -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fno-strict-aliasing -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/include/python2.7 -c tdlib.cpp -o tdlib.o

pylib:
	$(CC) -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-z,relro -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -Wl,-z,relro -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security tdlib.o -L. -ltd -o tdlib.so

partial_clean:
	-@rm -rf *.o *.pyc 2>/dev/null || true
	-@rm tdlib.cpp libtd.so 2>/dev/null || true

test:
	python test.py

clean: 
	-@rm -rf *.o *.so  *.pyc 2>/dev/null || true
	-@rm tdlib.cpp libtd.so 2>/dev/null || true 
