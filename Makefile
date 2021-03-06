#############################################################################
# Makefile for building: Act_model
# Generated by qmake (2.01a) (Qt 4.8.7) on: do mrt. 5 14:18:26 2020
# Project:  CellularPotts2.pro
# Template: app
# Command: /usr/bin/qmake-qt4 QT\ +=\ qt3support -o Makefile CellularPotts2.pro
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DQT_QT3SUPPORT_LIB -DQT3_SUPPORT -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
CFLAGS        = -m64 -pipe -g -Wall -W -D_REENTRANT $(DEFINES)
CXXFLAGS      = -m64 -pipe -g -DQTGRAPHICS -Wall -W -D_REENTRANT $(DEFINES)
INCPATH       = -I/usr/share/qt4/mkspecs/linux-g++-64 -I. -I/usr/include/qt4/QtCore -I/usr/include/qt4/QtGui -I/usr/include/qt4/Qt3Support -I/usr/include/qt4 -I.
LINK          = g++
LFLAGS        = -no-pie
LIBS          = $(SUBLIBS)  -L/usr/lib/x86_64-linux-gnu -lQt3Support -lQtGui -lQtCore -lpthread 
AR            = ar cqs
RANLIB        = 
QMAKE         = /usr/bin/qmake-qt4
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = ca.cpp \
		hull.cpp \
		cell.cpp \
		conrec.cpp \
		dish.cpp \
		info.cpp \
		misc.cpp \
		output.cpp \
		parameter.cpp \
		parse.cpp \
		pde.cpp \
		random.cpp \
		crash.cpp \
		warning.cpp \
		Act_model.cpp \
		qtgraph.cpp moc_qtgraph.cpp
OBJECTS       = ca.o \
		hull.o \
		cell.o \
		conrec.o \
		dish.o \
		info.o \
		misc.o \
		output.o \
		parameter.o \
		parse.o \
		pde.o \
		random.o \
		crash.o \
		warning.o \
		Act_model.o \
		qtgraph.o \
		moc_qtgraph.o
DIST          = /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/debug.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/shared.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		CellularPotts2.pro
QMAKE_TARGET  = Act_model
DESTDIR       = 
TARGET        = Act_model

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)
	{ test -n "$(DESTDIR)" && DESTDIR="$(DESTDIR)" || DESTDIR=.; } && test $$(gdb --version | sed -e 's,[^0-9][^0-9]*\([0-9]\)\.\([0-9]\).*,\1\2,;q') -gt 72 && gdb --nx --batch --quiet -ex 'set confirm off' -ex "save gdb-index $$DESTDIR" -ex quit '$(TARGET)' && test -f $(TARGET).gdb-index && objcopy --add-section '.gdb_index=$(TARGET).gdb-index' --set-section-flags '.gdb_index=readonly' '$(TARGET)' '$(TARGET)' && rm -f $(TARGET).gdb-index || true

Makefile: CellularPotts2.pro  /usr/share/qt4/mkspecs/linux-g++-64/qmake.conf /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/debug.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/shared.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		/usr/lib/x86_64-linux-gnu/libQt3Support.prl \
		/usr/lib/x86_64-linux-gnu/libQtGui.prl \
		/usr/lib/x86_64-linux-gnu/libQtCore.prl
	$(QMAKE) QT\ +=\ qt3support -o Makefile CellularPotts2.pro
/usr/share/qt4/mkspecs/common/unix.conf:
/usr/share/qt4/mkspecs/common/linux.conf:
/usr/share/qt4/mkspecs/common/gcc-base.conf:
/usr/share/qt4/mkspecs/common/gcc-base-unix.conf:
/usr/share/qt4/mkspecs/common/g++-base.conf:
/usr/share/qt4/mkspecs/common/g++-unix.conf:
/usr/share/qt4/mkspecs/qconfig.pri:
/usr/share/qt4/mkspecs/features/qt_functions.prf:
/usr/share/qt4/mkspecs/features/qt_config.prf:
/usr/share/qt4/mkspecs/features/exclusive_builds.prf:
/usr/share/qt4/mkspecs/features/default_pre.prf:
/usr/share/qt4/mkspecs/features/debug.prf:
/usr/share/qt4/mkspecs/features/default_post.prf:
/usr/share/qt4/mkspecs/features/shared.prf:
/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf:
/usr/share/qt4/mkspecs/features/warn_on.prf:
/usr/share/qt4/mkspecs/features/qt.prf:
/usr/share/qt4/mkspecs/features/unix/thread.prf:
/usr/share/qt4/mkspecs/features/moc.prf:
/usr/share/qt4/mkspecs/features/resources.prf:
/usr/share/qt4/mkspecs/features/uic.prf:
/usr/share/qt4/mkspecs/features/yacc.prf:
/usr/share/qt4/mkspecs/features/lex.prf:
/usr/share/qt4/mkspecs/features/include_source_dir.prf:
/usr/lib/x86_64-linux-gnu/libQt3Support.prl:
/usr/lib/x86_64-linux-gnu/libQtGui.prl:
/usr/lib/x86_64-linux-gnu/libQtCore.prl:
qmake:  FORCE
	@$(QMAKE) QT\ +=\ qt3support -o Makefile CellularPotts2.pro

dist: 
	@$(CHK_DIR_EXISTS) .tmp/Act_model1.0.0 || $(MKDIR) .tmp/Act_model1.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) .tmp/Act_model1.0.0/ && $(COPY_FILE) --parents ca.h hull.h cell.h conrec.h dish.h graph.h info.h misc.h output.h parameter.h parse.h pde.h random.h sqr.h sticky.h crash.h warning.h qtgraph.h .tmp/Act_model1.0.0/ && $(COPY_FILE) --parents ca.cpp hull.cpp cell.cpp conrec.cpp dish.cpp info.cpp misc.cpp output.cpp parameter.cpp parse.cpp pde.cpp random.cpp crash.cpp warning.cpp Act_model.cpp qtgraph.cpp .tmp/Act_model1.0.0/ && (cd `dirname .tmp/Act_model1.0.0` && $(TAR) Act_model1.0.0.tar Act_model1.0.0 && $(COMPRESS) Act_model1.0.0.tar) && $(MOVE) `dirname .tmp/Act_model1.0.0`/Act_model1.0.0.tar.gz . && $(DEL_FILE) -r .tmp/Act_model1.0.0


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


check: first

mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compiler_moc_header_make_all: moc_qtgraph.cpp
compiler_moc_header_clean:
	-$(DEL_FILE) moc_qtgraph.cpp
moc_qtgraph.cpp: graph.h \
		qtgraph.h
	/usr/lib/x86_64-linux-gnu/qt4/bin/moc $(DEFINES) $(INCPATH) qtgraph.h -o moc_qtgraph.cpp

compiler_rcc_make_all:
compiler_rcc_clean:
compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_uic_make_all:
compiler_uic_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: compiler_moc_header_clean 

####### Compile

ca.o: ca.cpp sticky.h \
		random.h \
		ca.h \
		mainpage.h \
		graph.h \
		pde.h \
		cell.h \
		parameter.h \
		dish.h \
		sqr.h \
		crash.h \
		hull.h \
		1.xpm
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o ca.o ca.cpp

hull.o: hull.cpp hull.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o hull.o hull.cpp

cell.o: cell.cpp cell.h \
		parameter.h \
		sticky.h \
		dish.h \
		graph.h \
		random.h \
		pde.h \
		ca.h \
		mainpage.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o cell.o cell.cpp

conrec.o: conrec.cpp graph.h \
		conrec.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o conrec.o conrec.cpp

dish.o: dish.cpp dish.h \
		graph.h \
		random.h \
		pde.h \
		cell.h \
		parameter.h \
		ca.h \
		mainpage.h \
		sticky.h \
		info.h \
		crash.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dish.o dish.cpp

info.o: info.cpp dish.h \
		graph.h \
		random.h \
		pde.h \
		cell.h \
		parameter.h \
		ca.h \
		mainpage.h \
		info.h \
		misc.h \
		parse.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o info.o info.cpp

misc.o: misc.cpp sticky.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o misc.o misc.cpp

output.o: output.cpp warning.h \
		parameter.h \
		output.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o output.o output.cpp

parameter.o: parameter.cpp parameter.h \
		output.h \
		parse.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parameter.o parameter.cpp

parse.o: parse.cpp warning.h \
		parse.h \
		output.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parse.o parse.cpp

pde.o: pde.cpp crash.h \
		parameter.h \
		ca.h \
		mainpage.h \
		graph.h \
		pde.h \
		cell.h \
		conrec.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pde.o pde.cpp

random.o: random.cpp random.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o random.o random.cpp

crash.o: crash.cpp sticky.h \
		crash.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o crash.o crash.cpp

warning.o: warning.cpp warning.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o warning.o warning.cpp

Act_model.o: Act_model.cpp dish.h \
		graph.h \
		random.h \
		pde.h \
		cell.h \
		parameter.h \
		ca.h \
		mainpage.h \
		info.h \
		sqr.h \
		leonie.cpp \
		qtgraph.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Act_model.o Act_model.cpp

qtgraph.o: qtgraph.cpp qtgraph.h \
		graph.h \
		parameter.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o qtgraph.o qtgraph.cpp

moc_qtgraph.o: moc_qtgraph.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_qtgraph.o moc_qtgraph.cpp

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

