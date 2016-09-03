CC=g++
CFLAGS=-O3 -march=native -std=c++1y
LDFLAGS=-lboost_graph -lboost_thread -lboost_system
PDFLATEX = pdflatex
BIBTEX = bibtex

docdatadir = ${pkgdatadir}
docdata_DATA = documentation.pdf

# fallback to Results... sometimes
TEX_ENV=TEXINPUTS=$(builddir)/src:$(srcdir)/src:$(srcdir)/Results:

documentation.pdf: documentation.tex  documentation.bib
	$(TEX_ENV) $(PDFLATEX) $<
	$(BIBTEX) $<
	$(TEX_ENV) $(PDFLATEX) $<

SUBDIRS=src .

table1: ./src/table1.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table1 ./src/table1.cpp

table2: ./src/table2.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table2 ./src/table2.cpp

table3: ./src/table3.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table3 ./src/table3.cpp

table4: ./src/table4.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table4 ./src/table4.cpp

table5: ./src/table5.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table5 ./src/table5.cpp

table6: ./src/table6.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table6 ./src/table6.cpp

table7: ./src/table7.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table7 ./src/table7.cpp

table8: ./src/table8.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table8 ./src/table8.cpp

table9: ./src/table9.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table9 ./src/table9.cpp

table10: ./src/table10.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table10 ./src/table10.cpp

table11: ./src/table11.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table11 ./src/table11.cpp

table12: ./src/table12.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table12 ./src/table12.cpp

table13: ./src/table13.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table13 ./src/table13.cpp

table14: ./src/table14.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table14 ./src/table14.cpp

table15: ./src/table15.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table15 ./src/table15.cpp

table16: ./src/table16.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table16 ./src/table16.cpp

table17: ./src/table17.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table17 ./src/table17.cpp

table18: ./src/table18.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table18 ./src/table18.cpp

table19: ./src/table19.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table19 ./src/table19.cpp

table20: ./src/table20.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table20 ./src/table20.cpp

table21: ./src/table21.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table21 ./src/table21.cpp

table22: ./src/table22.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table22 ./src/table22.cpp

table23: ./src/table23.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table23 ./src/table23.cpp

table24: ./src/table24.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table24 ./src/table24.cpp

table25: ./src/table25.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table25 ./src/table25.cpp

table26: ./src/table26.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table26 ./src/table26.cpp

table27: ./src/table27.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table27 ./src/table27.cpp

table28: ./src/table28.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table28 ./src/table28.cpp

table29: ./src/table29.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table29 ./src/table29.cpp

table30: ./src/table30.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table30 ./src/table30.cpp

table31: ./src/table31.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table31 ./src/table31.cpp

table32: ./src/table32.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table32 ./src/table32.cpp

table33: ./src/table33.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table33 ./src/table33.cpp

table34: ./src/table34.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -o table34 ./src/table34.cpp

GRAPHDIR=$(srcdir)/../graphs

1_PP.tex: table1
	: > $@
	./table1 stdlib $(GRAPHDIR)/stdlib/*; \
	./table1 contiki $(GRAPHDIR)/contiki/*; \
	./table1 fuzix $(GRAPHDIR)/fuzix/*; \
	./table1 coremark $(GRAPHDIR)/coremark/*; \
	./table1 whetstone $(GRAPHDIR)/whetstone/*; \
	./table1 dhrystone $(GRAPHDIR)/dhrystone/*; \
	./table1 dimacs $(GRAPHDIR)/dimacs/*; \
	./table1 zoo $(GRAPHDIR)/zoo/*;

data2:  table2
	./table2 contiki $(GRAPHDIR)/contiki/*; \
	./table2 fuzix $(GRAPHDIR)/fuzix/*;

data3:  table3
	./table3 contiki $(GRAPHDIR)/contiki/*; \
	./table3 fuzix $(GRAPHDIR)/fuzix/*;

data4:  table4
	./table4 contiki $(GRAPHDIR)/contiki/*; \
	./table4 fuzix $(GRAPHDIR)/fuzix/*;

data5:  table5
	./table5 contiki $(GRAPHDIR)/contiki/*; \
	./table5 fuzix $(GRAPHDIR)/fuzix/*;

data6:  table6
	./table6 contiki $(GRAPHDIR)/contiki/*; \
	./table6 fuzix $(GRAPHDIR)/fuzix/*;

data9:  table9
	./table9 stdlib $(GRAPHDIR)/stdlib/*; \
	./table9 contiki $(GRAPHDIR)/contiki/*; \
	./table9 fuzix $(GRAPHDIR)/fuzix/*; \
	./table9 coremark $(GRAPHDIR)/coremark/*; \
	./table9 whetstone $(GRAPHDIR)/whetstone/*; \
	./table9 dhrystone $(GRAPHDIR)/dhrystone/*; \
	./table9 dimacs $(GRAPHDIR)/dimacs/*; \
	./table9 zoo $(GRAPHDIR)/zoo/*;

data10:  table10
	./table10 stdlib $(GRAPHDIR)/stdlib/*; \
	./table10 contiki $(GRAPHDIR)/contiki/*; \
	./table10 fuzix $(GRAPHDIR)/fuzix/*; \
	./table10 coremark $(GRAPHDIR)/coremark/*; \
	./table10 whetstone $(GRAPHDIR)/whetstone/*; \
	./table10 dhrystone $(GRAPHDIR)/dhrystone/*; \
	./table10 dimacs $(GRAPHDIR)/dimacs/*; \
	./table10 zoo $(GRAPHDIR)/zoo/*;

data11:  table11
	./table11 stdlib $(GRAPHDIR)/stdlib/*; \
	./table11 contiki $(GRAPHDIR)/contiki/*; \
	./table11 fuzix $(GRAPHDIR)/fuzix/*; \
	./table11 coremark $(GRAPHDIR)/coremark/*; \
	./table11 whetstone $(GRAPHDIR)/whetstone/*; \
	./table11 dhrystone $(GRAPHDIR)/dhrystone/*; \
	./table11 dimacs $(GRAPHDIR)/dimacs/*; \
	./table11 zoo $(GRAPHDIR)/zoo/*;

data12:  table12
	./table12 stdlib $(GRAPHDIR)/stdlib/*; \
	./table12 contiki $(GRAPHDIR)/contiki/*; \
	./table12 fuzix $(GRAPHDIR)/fuzix/*; \
	./table12 coremark $(GRAPHDIR)/coremark/*; \
	./table12 whetstone $(GRAPHDIR)/whetstone/*; \
	./table12 dhrystone $(GRAPHDIR)/dhrystone/*; \
	./table12 dimacs $(GRAPHDIR)/dimacs/*; \
	./table12 zoo $(GRAPHDIR)/zoo/*;

data13:  table13
	./table13 stdlib $(GRAPHDIR)/stdlib/*; \
	./table13 contiki $(GRAPHDIR)/contiki/*; \
	./table13 fuzix $(GRAPHDIR)/fuzix/*; \
	./table13 coremark $(GRAPHDIR)/coremark/*; \
	./table13 whetstone $(GRAPHDIR)/whetstone/*; \
	./table13 dhrystone $(GRAPHDIR)/dhrystone/*;

data14:  table14
	./table14 stdlib $(GRAPHDIR)/stdlib/*; \
	./table14 contiki $(GRAPHDIR)/contiki/*; \
	./table14 fuzix $(GRAPHDIR)/fuzix/*; \
	./table14 coremark $(GRAPHDIR)/coremark/*; \
	./table14 whetstone $(GRAPHDIR)/whetstone/*; \
	./table14 dhrystone $(GRAPHDIR)/dhrystone/*;

data15:  table15
	./table15 stdlib $(GRAPHDIR)/stdlib/*; \
	./table15 contiki $(GRAPHDIR)/contiki/*; \
	./table15 fuzix $(GRAPHDIR)/fuzix/*; \
	./table15 coremark $(GRAPHDIR)/coremark/*; \
	./table15 whetstone $(GRAPHDIR)/whetstone/*; \
	./table15 dhrystone $(GRAPHDIR)/dhrystone/*;

data16:  table16
	./table16 stdlib $(GRAPHDIR)/stdlib/*; \
	./table16 contiki $(GRAPHDIR)/contiki/*; \
	./table16 fuzix $(GRAPHDIR)/fuzix/*; \
	./table16 coremark $(GRAPHDIR)/coremark/*; \
	./table16 whetstone $(GRAPHDIR)/whetstone/*; \
	./table16 dhrystone $(GRAPHDIR)/dhrystone/*;

data17:  table17
	./table17 stdlib $(GRAPHDIR)/stdlib/*; \
	./table17 contiki $(GRAPHDIR)/contiki/*; \
	./table17 fuzix $(GRAPHDIR)/fuzix/*; \
	./table17 coremark $(GRAPHDIR)/coremark/*; \
	./table17 whetstone $(GRAPHDIR)/whetstone/*; \
	./table17 dhrystone $(GRAPHDIR)/dhrystone/*;

data18:  table18
	./table18 stdlib $(GRAPHDIR)/stdlib/*; \
	./table18 contiki $(GRAPHDIR)/contiki/*; \
	./table18 fuzix $(GRAPHDIR)/fuzix/*; \
	./table18 coremark $(GRAPHDIR)/coremark/*; \
	./table18 whetstone $(GRAPHDIR)/whetstone/*; \
	./table18 dhrystone $(GRAPHDIR)/dhrystone/*;

data19:  table19
	./table19 stdlib $(GRAPHDIR)/stdlib/*; \
	./table19 contiki $(GRAPHDIR)/contiki/*; \
	./table19 fuzix $(GRAPHDIR)/fuzix/*; \
	./table19 coremark $(GRAPHDIR)/coremark/*; \
	./table19 whetstone $(GRAPHDIR)/whetstone/*; \
	./table19 dhrystone $(GRAPHDIR)/dhrystone/*;

data20:  table20
	./table20 stdlib $(GRAPHDIR)/stdlib/*; \
	./table20 contiki $(GRAPHDIR)/contiki/*; \
	./table20 fuzix $(GRAPHDIR)/fuzix/*; \
	./table20 coremark $(GRAPHDIR)/coremark/*; \
	./table20 whetstone $(GRAPHDIR)/whetstone/*; \
	./table20 dhrystone $(GRAPHDIR)/dhrystone/*;

data23:  table23
	./table23 stdlib $(GRAPHDIR)/stdlib/*; \
	./table23 contiki $(GRAPHDIR)/contiki/*; \
	./table23 fuzix $(GRAPHDIR)/fuzix/*; \
	./table23 coremark $(GRAPHDIR)/coremark/*; \
	./table23 whetstone $(GRAPHDIR)/whetstone/*; \
	./table23 dhrystone $(GRAPHDIR)/dhrystone/*;

data24:  table24
	./table24 stdlib $(GRAPHDIR)/stdlib/*; \
	./table24 contiki $(GRAPHDIR)/contiki/*; \
	./table24 fuzix $(GRAPHDIR)/fuzix/*; \
	./table24 coremark $(GRAPHDIR)/coremark/*; \
	./table24 whetstone $(GRAPHDIR)/whetstone/*; \
	./table24 dhrystone $(GRAPHDIR)/dhrystone/*;

data25:  table25
	./table25 stdlib $(GRAPHDIR)/stdlib/*; \
	./table25 contiki $(GRAPHDIR)/contiki/*; \
	./table25 fuzix $(GRAPHDIR)/fuzix/*; \
	./table25 coremark $(GRAPHDIR)/coremark/*; \
	./table25 whetstone $(GRAPHDIR)/whetstone/*; \
	./table25 dhrystone $(GRAPHDIR)/dhrystone/*;

data26:  table26
	./table26 stdlib $(GRAPHDIR)/stdlib/*; \
	./table26 contiki $(GRAPHDIR)/contiki/*; \
	./table26 fuzix $(GRAPHDIR)/fuzix/*; \
	./table26 coremark $(GRAPHDIR)/coremark/*; \
	./table26 whetstone $(GRAPHDIR)/whetstone/*; \
	./table26 dhrystone $(GRAPHDIR)/dhrystone/*;

data27:  table27
	./table27 stdlib $(GRAPHDIR)/stdlib/*; \
	./table27 contiki $(GRAPHDIR)/contiki/*; \
	./table27 fuzix $(GRAPHDIR)/fuzix/*; \
	./table27 coremark $(GRAPHDIR)/coremark/*; \
	./table27 whetstone $(GRAPHDIR)/whetstone/*; \
	./table27 dhrystone $(GRAPHDIR)/dhrystone/*;

data28:  table28
	./table28 stdlib $(GRAPHDIR)/stdlib/*; \
	./table28 contiki $(GRAPHDIR)/contiki/*; \
	./table28 fuzix $(GRAPHDIR)/fuzix/*; \
	./table28 coremark $(GRAPHDIR)/coremark/*; \
	./table28 whetstone $(GRAPHDIR)/whetstone/*; \
	./table28 dhrystone $(GRAPHDIR)/dhrystone/*;

data29:  table29
	./table29 stdlib $(GRAPHDIR)/stdlib/*; \
	./table29 contiki $(GRAPHDIR)/contiki/*; \
	./table29 fuzix $(GRAPHDIR)/fuzix/*; \
	./table29 coremark $(GRAPHDIR)/coremark/*; \
	./table29 whetstone $(GRAPHDIR)/whetstone/*; \
	./table29 dhrystone $(GRAPHDIR)/dhrystone/*;

data30:  table30
	./table30 stdlib $(GRAPHDIR)/stdlib/*; \
	./table30 contiki $(GRAPHDIR)/contiki/*; \
	./table30 fuzix $(GRAPHDIR)/fuzix/*; \
	./table30 coremark $(GRAPHDIR)/coremark/*; \
	./table30 whetstone $(GRAPHDIR)/whetstone/*; \
	./table30 dhrystone $(GRAPHDIR)/dhrystone/*;

data31:  table31
	./table31 stdlib $(GRAPHDIR)/stdlib/*; \
	./table31 contiki $(GRAPHDIR)/contiki/*; \
	./table31 fuzix $(GRAPHDIR)/fuzix/*; \
	./table31 coremark $(GRAPHDIR)/coremark/*; \
	./table31 whetstone $(GRAPHDIR)/whetstone/*; \
	./table31 dhrystone $(GRAPHDIR)/dhrystone/*;

data32:  table32
	./table32 stdlib $(GRAPHDIR)/stdlib/*; \
	./table32 contiki $(GRAPHDIR)/contiki/*; \
	./table32 fuzix $(GRAPHDIR)/fuzix/*; \
	./table32 coremark $(GRAPHDIR)/coremark/*; \
	./table32 whetstone $(GRAPHDIR)/whetstone/*; \
	./table32 dhrystone $(GRAPHDIR)/dhrystone/*;

data33:  table33
	./table33 stdlib $(GRAPHDIR)/stdlib/*; \
	./table33 contiki $(GRAPHDIR)/contiki/*; \
	./table33 fuzix $(GRAPHDIR)/fuzix/*; \
	./table33 coremark $(GRAPHDIR)/coremark/*; \
	./table33 whetstone $(GRAPHDIR)/whetstone/*; \
	./table33 dhrystone $(GRAPHDIR)/dhrystone/*;

data34:  table34
	./table34 stdlib $(GRAPHDIR)/stdlib/*; \
	./table34 contiki $(GRAPHDIR)/contiki/*; \
	./table34 fuzix $(GRAPHDIR)/fuzix/*; \
	./table34 coremark $(GRAPHDIR)/coremark/*; \
	./table34 whetstone $(GRAPHDIR)/whetstone/*; \
	./table34 dhrystone $(GRAPHDIR)/dhrystone/*;

CLEANFILES=*.aux *.blg *.log *.out *.pdf *.toc table*
