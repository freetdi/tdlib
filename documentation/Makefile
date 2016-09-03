CC=g++
CFLAGS=-O3 -march=native -std=c++1y
LDFLAGS=-lboost_graph -lboost_thread -lboost_system

documentation:	documentation.tex  documentation.bib
	pdflatex documentation.tex; bibtex documentation.aux; pdflatex documentation.tex

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

data1:  table1
	./table1 stdlib ./graphs/stdlib/*; \
	./table1 contiki ./graphs/contiki/*; \
	./table1 fuzix ./graphs/fuzix/*; \
	./table1 coremark ./graphs/coremark/*; \
	./table1 whetstone ./graphs/whetstone/*; \
	./table1 dhrystone ./graphs/dhrystone/*; \
	./table1 dimacs ./graphs/dimacs/*; \
	./table1 zoo ./graphs/zoo/*;

data2:  table2
	./table2 contiki ./graphs/contiki/*; \
	./table2 fuzix ./graphs/fuzix/*;

data3:  table3
	./table3 contiki ./graphs/contiki/*; \
	./table3 fuzix ./graphs/fuzix/*;

data4:  table4
	./table4 contiki ./graphs/contiki/*; \
	./table4 fuzix ./graphs/fuzix/*;

data5:  table5
	./table5 contiki ./graphs/contiki/*; \
	./table5 fuzix ./graphs/fuzix/*;

data6:  table6
	./table6 contiki ./graphs/contiki/*; \
	./table6 fuzix ./graphs/fuzix/*;

data9:  table9
	./table9 stdlib ./graphs/stdlib/*; \
	./table9 contiki ./graphs/contiki/*; \
	./table9 fuzix ./graphs/fuzix/*; \
	./table9 coremark ./graphs/coremark/*; \
	./table9 whetstone ./graphs/whetstone/*; \
	./table9 dhrystone ./graphs/dhrystone/*; \
	./table9 dimacs ./graphs/dimacs/*; \
	./table9 zoo ./graphs/zoo/*;

data10:  table10
	./table10 stdlib ./graphs/stdlib/*; \
	./table10 contiki ./graphs/contiki/*; \
	./table10 fuzix ./graphs/fuzix/*; \
	./table10 coremark ./graphs/coremark/*; \
	./table10 whetstone ./graphs/whetstone/*; \
	./table10 dhrystone ./graphs/dhrystone/*; \
	./table10 dimacs ./graphs/dimacs/*; \
	./table10 zoo ./graphs/zoo/*;

data11:  table11
	./table11 stdlib ./graphs/stdlib/*; \
	./table11 contiki ./graphs/contiki/*; \
	./table11 fuzix ./graphs/fuzix/*; \
	./table11 coremark ./graphs/coremark/*; \
	./table11 whetstone ./graphs/whetstone/*; \
	./table11 dhrystone ./graphs/dhrystone/*; \
	./table11 dimacs ./graphs/dimacs/*; \
	./table11 zoo ./graphs/zoo/*;

data12:  table12
	./table12 stdlib ./graphs/stdlib/*; \
	./table12 contiki ./graphs/contiki/*; \
	./table12 fuzix ./graphs/fuzix/*; \
	./table12 coremark ./graphs/coremark/*; \
	./table12 whetstone ./graphs/whetstone/*; \
	./table12 dhrystone ./graphs/dhrystone/*; \
	./table12 dimacs ./graphs/dimacs/*; \
	./table12 zoo ./graphs/zoo/*;

data13:  table13
	./table13 stdlib ./graphs/stdlib/*; \
	./table13 contiki ./graphs/contiki/*; \
	./table13 fuzix ./graphs/fuzix/*; \
	./table13 coremark ./graphs/coremark/*; \
	./table13 whetstone ./graphs/whetstone/*; \
	./table13 dhrystone ./graphs/dhrystone/*;

data14:  table14
	./table14 stdlib ./graphs/stdlib/*; \
	./table14 contiki ./graphs/contiki/*; \
	./table14 fuzix ./graphs/fuzix/*; \
	./table14 coremark ./graphs/coremark/*; \
	./table14 whetstone ./graphs/whetstone/*; \
	./table14 dhrystone ./graphs/dhrystone/*;

data15:  table15
	./table15 stdlib ./graphs/stdlib/*; \
	./table15 contiki ./graphs/contiki/*; \
	./table15 fuzix ./graphs/fuzix/*; \
	./table15 coremark ./graphs/coremark/*; \
	./table15 whetstone ./graphs/whetstone/*; \
	./table15 dhrystone ./graphs/dhrystone/*;

data16:  table16
	./table16 stdlib ./graphs/stdlib/*; \
	./table16 contiki ./graphs/contiki/*; \
	./table16 fuzix ./graphs/fuzix/*; \
	./table16 coremark ./graphs/coremark/*; \
	./table16 whetstone ./graphs/whetstone/*; \
	./table16 dhrystone ./graphs/dhrystone/*;

data17:  table17
	./table17 stdlib ./graphs/stdlib/*; \
	./table17 contiki ./graphs/contiki/*; \
	./table17 fuzix ./graphs/fuzix/*; \
	./table17 coremark ./graphs/coremark/*; \
	./table17 whetstone ./graphs/whetstone/*; \
	./table17 dhrystone ./graphs/dhrystone/*;

data18:  table18
	./table18 stdlib ./graphs/stdlib/*; \
	./table18 contiki ./graphs/contiki/*; \
	./table18 fuzix ./graphs/fuzix/*; \
	./table18 coremark ./graphs/coremark/*; \
	./table18 whetstone ./graphs/whetstone/*; \
	./table18 dhrystone ./graphs/dhrystone/*;

data19:  table19
	./table19 stdlib ./graphs/stdlib/*; \
	./table19 contiki ./graphs/contiki/*; \
	./table19 fuzix ./graphs/fuzix/*; \
	./table19 coremark ./graphs/coremark/*; \
	./table19 whetstone ./graphs/whetstone/*; \
	./table19 dhrystone ./graphs/dhrystone/*;

data20:  table20
	./table20 stdlib ./graphs/stdlib/*; \
	./table20 contiki ./graphs/contiki/*; \
	./table20 fuzix ./graphs/fuzix/*; \
	./table20 coremark ./graphs/coremark/*; \
	./table20 whetstone ./graphs/whetstone/*; \
	./table20 dhrystone ./graphs/dhrystone/*;

data23:  table23
	./table23 stdlib ./graphs/stdlib/*; \
	./table23 contiki ./graphs/contiki/*; \
	./table23 fuzix ./graphs/fuzix/*; \
	./table23 coremark ./graphs/coremark/*; \
	./table23 whetstone ./graphs/whetstone/*; \
	./table23 dhrystone ./graphs/dhrystone/*;

data24:  table24
	./table24 stdlib ./graphs/stdlib/*; \
	./table24 contiki ./graphs/contiki/*; \
	./table24 fuzix ./graphs/fuzix/*; \
	./table24 coremark ./graphs/coremark/*; \
	./table24 whetstone ./graphs/whetstone/*; \
	./table24 dhrystone ./graphs/dhrystone/*;

data25:  table25
	./table25 stdlib ./graphs/stdlib/*; \
	./table25 contiki ./graphs/contiki/*; \
	./table25 fuzix ./graphs/fuzix/*; \
	./table25 coremark ./graphs/coremark/*; \
	./table25 whetstone ./graphs/whetstone/*; \
	./table25 dhrystone ./graphs/dhrystone/*;

data26:  table26
	./table26 stdlib ./graphs/stdlib/*; \
	./table26 contiki ./graphs/contiki/*; \
	./table26 fuzix ./graphs/fuzix/*; \
	./table26 coremark ./graphs/coremark/*; \
	./table26 whetstone ./graphs/whetstone/*; \
	./table26 dhrystone ./graphs/dhrystone/*;

data27:  table27
	./table27 stdlib ./graphs/stdlib/*; \
	./table27 contiki ./graphs/contiki/*; \
	./table27 fuzix ./graphs/fuzix/*; \
	./table27 coremark ./graphs/coremark/*; \
	./table27 whetstone ./graphs/whetstone/*; \
	./table27 dhrystone ./graphs/dhrystone/*;

data28:  table28
	./table28 stdlib ./graphs/stdlib/*; \
	./table28 contiki ./graphs/contiki/*; \
	./table28 fuzix ./graphs/fuzix/*; \
	./table28 coremark ./graphs/coremark/*; \
	./table28 whetstone ./graphs/whetstone/*; \
	./table28 dhrystone ./graphs/dhrystone/*;

data29:  table29
	./table29 stdlib ./graphs/stdlib/*; \
	./table29 contiki ./graphs/contiki/*; \
	./table29 fuzix ./graphs/fuzix/*; \
	./table29 coremark ./graphs/coremark/*; \
	./table29 whetstone ./graphs/whetstone/*; \
	./table29 dhrystone ./graphs/dhrystone/*;

data30:  table30
	./table30 stdlib ./graphs/stdlib/*; \
	./table30 contiki ./graphs/contiki/*; \
	./table30 fuzix ./graphs/fuzix/*; \
	./table30 coremark ./graphs/coremark/*; \
	./table30 whetstone ./graphs/whetstone/*; \
	./table30 dhrystone ./graphs/dhrystone/*;

data31:  table31
	./table31 stdlib ./graphs/stdlib/*; \
	./table31 contiki ./graphs/contiki/*; \
	./table31 fuzix ./graphs/fuzix/*; \
	./table31 coremark ./graphs/coremark/*; \
	./table31 whetstone ./graphs/whetstone/*; \
	./table31 dhrystone ./graphs/dhrystone/*;

data32:  table32
	./table32 stdlib ./graphs/stdlib/*; \
	./table32 contiki ./graphs/contiki/*; \
	./table32 fuzix ./graphs/fuzix/*; \
	./table32 coremark ./graphs/coremark/*; \
	./table32 whetstone ./graphs/whetstone/*; \
	./table32 dhrystone ./graphs/dhrystone/*;

data33:  table33
	./table33 stdlib ./graphs/stdlib/*; \
	./table33 contiki ./graphs/contiki/*; \
	./table33 fuzix ./graphs/fuzix/*; \
	./table33 coremark ./graphs/coremark/*; \
	./table33 whetstone ./graphs/whetstone/*; \
	./table33 dhrystone ./graphs/dhrystone/*;

data34:  table34
	./table34 stdlib ./graphs/stdlib/*; \
	./table34 contiki ./graphs/contiki/*; \
	./table34 fuzix ./graphs/fuzix/*; \
	./table34 coremark ./graphs/coremark/*; \
	./table34 whetstone ./graphs/whetstone/*; \
	./table34 dhrystone ./graphs/dhrystone/*;


clean:
	rm -f documentation.aux documentation.blg documentation.log documentation.out documentation.pdf documentation.toc table* ./Results/*
