AM_DEFAULT_SOURCE_EXT = .cpp

BUILT_SOURCES = tdlib

#hack: there's no <tdlib ...>, use -I. and this.
tdlib:
	$(LN_S) $(abs_top_srcdir)/src $@

LDADD = -lboost_graph -lboost_thread -lboost_system
CPPFLAGS = -I. -DRESULTS_HERE

noinst_PROGRAMS = \
	table1 \
	table2 \
	table3 \
	table4 \
	table5 \
	table6 \
	table7 \
	table8 \
	table9 \
	table10 \
	table11 \
	table12 \
	table13 \
	table14 \
	table15 \
	table16 \
	table17 \
	table18 \
	table19 \
	table20 \
	table21 \
	table22 \
	table23 \
	table24 \
	table25 \
	table26 \
	table27 \
	table28 \
	table29 \
	table30 \
	table31 \
	table32 \
	table33 \
	table34

BUILT_SOURCES += \
	1_PP.tex \
	2_lb1.tex \
	3_lb2.tex \
	4_lb_improved1.tex \
	5_lb_improved2.tex \
	6_lb_mcs.tex \
	9_MD.tex \
	10_PP_MD.tex \
	11_FI.tex \
	12_PP_FI.tex \
	13_thorup.tex \
	14_separator.tex \
	15_msvs_trivial.tex \
	16_sep_msvs.tex \
	17_FI_TM.tex \
	18_PP_FI_TM.tex \
	19_sep_TM.tex \
	20_separator.tex \
	23_max_clique.tex \
	24_max_independent_set.tex \
	25_min_vertex_cover.tex \
	26_min_dominating_set.tex \
	27_min_coloring.tex \
	28_sep_msvs.tex \
	29_sep_TM.tex \
	30_MD.tex \
	31_PP_MD.tex \
	32_FI.tex \
	33_FI_TM.tex \
	34_PP_FI_TM.tex

GRAPHDIR=$(srcdir)/../graphs

1_PP.tex: table1
	: > $@
	./table1 stdlib $(GRAPHDIR)/stdlib/*
	./table1 contiki $(GRAPHDIR)/contiki/*
	./table1 fuzix $(GRAPHDIR)/fuzix/*
	./table1 coremark $(GRAPHDIR)/coremark/*
	./table1 whetstone $(GRAPHDIR)/whetstone/*
	./table1 dhrystone $(GRAPHDIR)/dhrystone/*
	./table1 dimacs $(GRAPHDIR)/dimacs/*
	./table1 zoo $(GRAPHDIR)/zoo/*

2_lb1.tex: table2
	./table2 contiki $(GRAPHDIR)/contiki/*
	./table2 fuzix $(GRAPHDIR)/fuzix/*

3_lb2.tex: table3
	: > $@
	./table3 contiki $(GRAPHDIR)/contiki/*
	./table3 fuzix $(GRAPHDIR)/fuzix/*

4_lb_improved1.tex: table4
	: > $@
	./table4 contiki $(GRAPHDIR)/contiki/*
	./table4 fuzix $(GRAPHDIR)/fuzix/*

5_lb_improved2.tex: table5
	: > $@
	./table5 contiki $(GRAPHDIR)/contiki/*
	./table5 fuzix $(GRAPHDIR)/fuzix/*

6_lb_mcs.tex: table6
	: > $@
	./table6 contiki $(GRAPHDIR)/contiki/*
	./table6 fuzix $(GRAPHDIR)/fuzix/*

9_MD.tex: table9
	: > $@
	./table9 stdlib $(GRAPHDIR)/stdlib/*
	./table9 contiki $(GRAPHDIR)/contiki/*
	./table9 fuzix $(GRAPHDIR)/fuzix/*
	./table9 coremark $(GRAPHDIR)/coremark/*
	./table9 whetstone $(GRAPHDIR)/whetstone/*
	./table9 dhrystone $(GRAPHDIR)/dhrystone/*
	./table9 dimacs $(GRAPHDIR)/dimacs/*
	./table9 zoo $(GRAPHDIR)/zoo/*

10_PP_MD.tex: table10
	: > $@
	./table10 stdlib $(GRAPHDIR)/stdlib/*
	./table10 contiki $(GRAPHDIR)/contiki/*
	./table10 fuzix $(GRAPHDIR)/fuzix/*
	./table10 coremark $(GRAPHDIR)/coremark/*
	./table10 whetstone $(GRAPHDIR)/whetstone/*
	./table10 dhrystone $(GRAPHDIR)/dhrystone/*
	./table10 dimacs $(GRAPHDIR)/dimacs/*
	./table10 zoo $(GRAPHDIR)/zoo/*

11_FI.tex: table11
	: > $@
	./table11 stdlib $(GRAPHDIR)/stdlib/*
	./table11 contiki $(GRAPHDIR)/contiki/*
	./table11 fuzix $(GRAPHDIR)/fuzix/*
	./table11 coremark $(GRAPHDIR)/coremark/*
	./table11 whetstone $(GRAPHDIR)/whetstone/*
	./table11 dhrystone $(GRAPHDIR)/dhrystone/*
	./table11 dimacs $(GRAPHDIR)/dimacs/*
	./table11 zoo $(GRAPHDIR)/zoo/*

12_PP_FI.tex: table12
	: > $@
	./table12 stdlib $(GRAPHDIR)/stdlib/*
	./table12 contiki $(GRAPHDIR)/contiki/*
	./table12 fuzix $(GRAPHDIR)/fuzix/*
	./table12 coremark $(GRAPHDIR)/coremark/*
	./table12 whetstone $(GRAPHDIR)/whetstone/*
	./table12 dhrystone $(GRAPHDIR)/dhrystone/*
	./table12 dimacs $(GRAPHDIR)/dimacs/*
	./table12 zoo $(GRAPHDIR)/zoo/*

13_thorup.tex: table13
	: > $@
	./table13 stdlib $(GRAPHDIR)/stdlib/*
	./table13 contiki $(GRAPHDIR)/contiki/*
	./table13 fuzix $(GRAPHDIR)/fuzix/*
	./table13 coremark $(GRAPHDIR)/coremark/*
	./table13 whetstone $(GRAPHDIR)/whetstone/*
	./table13 dhrystone $(GRAPHDIR)/dhrystone/*

14_separator.tex: table14
	: > $@
	./table14 stdlib $(GRAPHDIR)/stdlib/*
	./table14 contiki $(GRAPHDIR)/contiki/*
	./table14 fuzix $(GRAPHDIR)/fuzix/*
	./table14 coremark $(GRAPHDIR)/coremark/*
	./table14 whetstone $(GRAPHDIR)/whetstone/*
	./table14 dhrystone $(GRAPHDIR)/dhrystone/*

15_msvs_trivial.tex: table15
	: > $@
	./table15 stdlib $(GRAPHDIR)/stdlib/*
	./table15 contiki $(GRAPHDIR)/contiki/*
	./table15 fuzix $(GRAPHDIR)/fuzix/*
	./table15 coremark $(GRAPHDIR)/coremark/*
	./table15 whetstone $(GRAPHDIR)/whetstone/*
	./table15 dhrystone $(GRAPHDIR)/dhrystone/*

16_sep_msvs.tex: table16
	: > $@
	./table16 stdlib $(GRAPHDIR)/stdlib/*
	./table16 contiki $(GRAPHDIR)/contiki/*
	./table16 fuzix $(GRAPHDIR)/fuzix/*
	./table16 coremark $(GRAPHDIR)/coremark/*
	./table16 whetstone $(GRAPHDIR)/whetstone/*
	./table16 dhrystone $(GRAPHDIR)/dhrystone/*

17_FI_TM.tex: table17
	: > $@
	./table17 stdlib $(GRAPHDIR)/stdlib/*
	./table17 contiki $(GRAPHDIR)/contiki/*
	./table17 fuzix $(GRAPHDIR)/fuzix/*
	./table17 coremark $(GRAPHDIR)/coremark/*
	./table17 whetstone $(GRAPHDIR)/whetstone/*
	./table17 dhrystone $(GRAPHDIR)/dhrystone/*

18_PP_FI_TM.tex: table18
	: > $@
	./table18 stdlib $(GRAPHDIR)/stdlib/*
	./table18 contiki $(GRAPHDIR)/contiki/*
	./table18 fuzix $(GRAPHDIR)/fuzix/*
	./table18 coremark $(GRAPHDIR)/coremark/*
	./table18 whetstone $(GRAPHDIR)/whetstone/*
	./table18 dhrystone $(GRAPHDIR)/dhrystone/*

19_sep_TM.tex: table19
	: > $@
	./table19 stdlib $(GRAPHDIR)/stdlib/*
	./table19 contiki $(GRAPHDIR)/contiki/*
	./table19 fuzix $(GRAPHDIR)/fuzix/*
	./table19 coremark $(GRAPHDIR)/coremark/*
	./table19 whetstone $(GRAPHDIR)/whetstone/*
	./table19 dhrystone $(GRAPHDIR)/dhrystone/*

20_separator.tex: table20
	: > $@
	./table20 stdlib $(GRAPHDIR)/stdlib/*
	./table20 contiki $(GRAPHDIR)/contiki/*
	./table20 fuzix $(GRAPHDIR)/fuzix/*
	./table20 coremark $(GRAPHDIR)/coremark/*
	./table20 whetstone $(GRAPHDIR)/whetstone/*
	./table20 dhrystone $(GRAPHDIR)/dhrystone/*

23_max_clique.tex: table23
	: > $@
	./table23 stdlib $(GRAPHDIR)/stdlib/*
	./table23 contiki $(GRAPHDIR)/contiki/*
	./table23 fuzix $(GRAPHDIR)/fuzix/*
	./table23 coremark $(GRAPHDIR)/coremark/*
	./table23 whetstone $(GRAPHDIR)/whetstone/*
	./table23 dhrystone $(GRAPHDIR)/dhrystone/*

24_max_independent_set.tex: table24
	: > $@
	./table24 stdlib $(GRAPHDIR)/stdlib/*
	./table24 contiki $(GRAPHDIR)/contiki/*
	./table24 fuzix $(GRAPHDIR)/fuzix/*
	./table24 coremark $(GRAPHDIR)/coremark/*
	./table24 whetstone $(GRAPHDIR)/whetstone/*
	./table24 dhrystone $(GRAPHDIR)/dhrystone/*

25_min_vertex_cover.tex: table25
	: > $@
	./table25 stdlib $(GRAPHDIR)/stdlib/*
	./table25 contiki $(GRAPHDIR)/contiki/*
	./table25 fuzix $(GRAPHDIR)/fuzix/*
	./table25 coremark $(GRAPHDIR)/coremark/*
	./table25 whetstone $(GRAPHDIR)/whetstone/*
	./table25 dhrystone $(GRAPHDIR)/dhrystone/*

26_min_dominating_set.tex: table26
	: > $@
	./table26 stdlib $(GRAPHDIR)/stdlib/*
	./table26 contiki $(GRAPHDIR)/contiki/*
	./table26 fuzix $(GRAPHDIR)/fuzix/*
	./table26 coremark $(GRAPHDIR)/coremark/*
	./table26 whetstone $(GRAPHDIR)/whetstone/*
	./table26 dhrystone $(GRAPHDIR)/dhrystone/*

27_min_coloring.tex: table27
	: > $@
	./table27 stdlib $(GRAPHDIR)/stdlib/*
	./table27 contiki $(GRAPHDIR)/contiki/*
	./table27 fuzix $(GRAPHDIR)/fuzix/*
	./table27 coremark $(GRAPHDIR)/coremark/*
	./table27 whetstone $(GRAPHDIR)/whetstone/*
	./table27 dhrystone $(GRAPHDIR)/dhrystone/*

28_sep_msvs.tex: table28
	: > $@
	./table28 stdlib $(GRAPHDIR)/stdlib/*
	./table28 contiki $(GRAPHDIR)/contiki/*
	./table28 fuzix $(GRAPHDIR)/fuzix/*
	./table28 coremark $(GRAPHDIR)/coremark/*
	./table28 whetstone $(GRAPHDIR)/whetstone/*
	./table28 dhrystone $(GRAPHDIR)/dhrystone/*

29_sep_TM.tex: table29
	: > $@
	./table29 stdlib $(GRAPHDIR)/stdlib/*
	./table29 contiki $(GRAPHDIR)/contiki/*
	./table29 fuzix $(GRAPHDIR)/fuzix/*
	./table29 coremark $(GRAPHDIR)/coremark/*
	./table29 whetstone $(GRAPHDIR)/whetstone/*
	./table29 dhrystone $(GRAPHDIR)/dhrystone/*

30_MD.tex: table30
	: > $@
	./table30 stdlib $(GRAPHDIR)/stdlib/*
	./table30 contiki $(GRAPHDIR)/contiki/*
	./table30 fuzix $(GRAPHDIR)/fuzix/*
	./table30 coremark $(GRAPHDIR)/coremark/*
	./table30 whetstone $(GRAPHDIR)/whetstone/*
	./table30 dhrystone $(GRAPHDIR)/dhrystone/*

31_PP_MD.tex: table31
	: > $@
	./table31 stdlib $(GRAPHDIR)/stdlib/*
	./table31 contiki $(GRAPHDIR)/contiki/*
	./table31 fuzix $(GRAPHDIR)/fuzix/*
	./table31 coremark $(GRAPHDIR)/coremark/*
	./table31 whetstone $(GRAPHDIR)/whetstone/*
	./table31 dhrystone $(GRAPHDIR)/dhrystone/*

32_FI.tex: table32
	: > $@
	./table32 stdlib $(GRAPHDIR)/stdlib/*
	./table32 contiki $(GRAPHDIR)/contiki/*
	./table32 fuzix $(GRAPHDIR)/fuzix/*
	./table32 coremark $(GRAPHDIR)/coremark/*
	./table32 whetstone $(GRAPHDIR)/whetstone/*
	./table32 dhrystone $(GRAPHDIR)/dhrystone/*

33_FI_TM.tex: table33
	: > $@
	./table33 stdlib $(GRAPHDIR)/stdlib/*
	./table33 contiki $(GRAPHDIR)/contiki/*
	./table33 fuzix $(GRAPHDIR)/fuzix/*
	./table33 coremark $(GRAPHDIR)/coremark/*
	./table33 whetstone $(GRAPHDIR)/whetstone/*
	./table33 dhrystone $(GRAPHDIR)/dhrystone/*

34_PP_FI_TM.tex: table34
	: > $@
	./table34 stdlib $(GRAPHDIR)/stdlib/*
	./table34 contiki $(GRAPHDIR)/contiki/*
	./table34 fuzix $(GRAPHDIR)/fuzix/*
	./table34 coremark $(GRAPHDIR)/coremark/*
	./table34 whetstone $(GRAPHDIR)/whetstone/*
	./table34 dhrystone $(GRAPHDIR)/dhrystone/*

EXTRA_DIST = tex.cpp helper.cpp
CLEANFILES=*.aux *.blg *.log *.out *.pdf *.toc table*
