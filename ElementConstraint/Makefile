#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: Makefile,v 1.1.1.1 2012/07/12 20:45:11 wku Exp $

#@file    Makefile
#@brief   Makefile for C++ example project using SCIP as a callable library
#@author  Thorsten Koch
#@author  Tobias Achterberg


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------
#VERSION        =        0.2a

SCIPDIR		=	lib/scip

############
# WKU: add CPoptimizer
#CONCERTDIR_CPO = /opt/ilog/cpx_studio_12.3/concert/
#CPLEXDIR_CPO = /opt/ilog/cpx_studio_12.3/cplex/

#ILOGINCFLAGS = -DIL_STD -O -DNDEBUG -I/opt/ilog/cpx_studio_12.3/cpoptimizer/include -I/opt/ilog/cpx_studio_12.3/concert/include -fPIC -fstrict-aliasing -pedantic -Wall -fexceptions -Wno-long-long -ffloat-store -m64 -DILOUSEMT -D_REENTRANT -DILM_REENTRANT

#ILOGLINKFLAGS = -L/opt/ilog/cpx_studio_12.3/cpoptimizer/lib/x86-64_sles10_4.1/static_pic -lcp -L/opt/ilog/cpx_studio_12.3/cplex/lib/x86-64_sles10_4.1/static_pic -lcplex -L/opt/ilog/cpx_studio_12.3/concert/lib/x86-64_sles10_4.1/static_pic -lconcert  -lpthread 

############
# WKU: ILOG scheduler
 #ILOGDIR		=	lib/ilog
 #CONCERTDIR=$(ILOGDIR)/concert29
 #SOLVERDIR=$(ILOGDIR)/solver67
 #SCHEDDIR=$(ILOGDIR)/scheduler67

 #CONCERTINCDIR = $(CONCERTDIR)/include
 #SOLVERINCDIR=$(SOLVERDIR)/include
 #SCHEDINCDIR=$(SCHEDDIR)/include
 
 #ILOGINCFLAGS=-fPIC -DIL_STD -I$(SCHEDINCDIR) -I$(SOLVERINCDIR) -I$(CONCERTINCDIR) 

 #SYSTEM     = x86-64_debian4.0_4.1
 #LIBFORMAT  = static_pic

 #CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
 #SOLVERLIBDIR  = $(SOLVERDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
 #SCHEDLIBDIR   = $(SCHEDDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

 #ILOGLINKFLAGS = -L$(SCHEDLIBDIR) -lschedule -L$(SOLVERLIBDIR) -lsolveriim -lsolver -lm -lpthread -L$(CONCERTLIBDIR) -lconcert
############

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	scip
MAINOBJ		=	cmain.o \
			cons_element.o \
			reader_ns.o \
			probdata_ns.o 						

MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.cppmain.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) lint/$(MAINNAME).lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: doc
doc:		
		cd doc; $(DOXY) $(MAINNAME).dxy

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):	
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:
ifneq ($(OBJDIR),)
		-rm -f $(OBJDIR)/*.o
		-rmdir -p $(OBJDIR)
endif
		-rm -f $(MAINFILE)
#.PHONY: test
#test:
#		cd check; \
#		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) #$(CONTINUE) $(LOCK) $(VERSION) $(LPS);

.PHONY: test
test:           $(MAINFILE)
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/getlastprob.awk);
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) "examples" $(LPS) $(VALGRIND);



.PHONY: lbbd_test
lbbd_test:
		cd check; \
		$(SHELL) ./checklbbd.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS);

.PHONY: tags
tags:
		rm -f TAGS; ctags -e -R -h ".c.cpp.h" --exclude=".*" src/ $(SCIPDIR)/src/; 

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.cpp|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
		$(ILOGLINKFLAGS) \
		$(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) \
		$(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX)\
		$(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(ILOGINCFLAGS) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
