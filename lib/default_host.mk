## $Id: Makefile 4384 2009-03-13 02:59:02Z yasuyuki $

#######################################################################
#                                                                     #
#                        Makefile ( GNU make )                        #
#                                                                     #
#######################################################################
TOP	?=

CPPSRC  ?=
CSRC    ?=
ASRC    ?=

#######################################################################
#                                                                     #
#         TOOL                                                        #
#                                                                     #
#######################################################################
CROSS	?=

AS       =$(CROSS)as
ASOPT   ?=

CC       =$(CROSS)clang++
#CCOPT   += -D__LINUX__

LD       =$(CROSS)clang++
LDOPT	?=


#######################################################################
#                                                                     #
#         PATH                                                        #
#                                                                     #
#######################################################################
INC_DIR += -I$(TOP)/.

OBJDIR ?= obj

#######################################################################
#                           Standard processing                       #
#######################################################################
CPPOBJS= ${patsubst %.cpp, ${OBJDIR}/%.o, ${notdir ${CPPSRC}}}
COBJS  = ${patsubst %.c, ${OBJDIR}/%.o, ${notdir ${CSRC}}}
AOBJS  = ${patsubst %.S, ${OBJDIR}/%.o, ${notdir ${ASRC}}}

OBJS   = ${CPPOBJS} ${COBJS} ${AOBJS}

all: make_kernel ${OBJDIR} ${OBJS} ${TARGET}

${TARGET}: ${OBJS}
	${CC} -o ${TARGET} ${OBJS} ${LIB_DIR} ${LIBS} ${LDOPT}

${OBJDIR}:
	@if [ ! -d ${OBJDIR} ]; then \
		echo ";; mkdir ${OBJDIR}"; mkdir ${OBJDIR}; \
	fi

$(CPPOBJS) : ${OBJDIR}/%.o : %.cpp
	${CC} ${CCOPT} ${INC_DIR} ${ADD_INC_DIR} -c $< -o $@

$(COBJS) : ${OBJDIR}/%.o : %.c
	${CC} ${CCOPT} ${INC_DIR} ${ADD_INC_DIR} -c $< -o $@

$(AOBJS) : ${OBJDIR}/%.o : %.S
	${CC} ${CCOPT} ${INC_DIR} ${ADD_INC_DIR} -c $< -o $@

.cpp.o:
	${CC} ${CCOPT} ${INC_DIR} -c $<	-o ${OBJDIR}/${notdir $@}

.c.o:
	${CC} ${CCOPT} ${INC_DIR} -c $<	-o ${OBJDIR}/${notdir $@}

.S.o:
	${CC} ${CCOPT} ${INC_DIR} -c $<	-o ${OBJDIR}/${notdir $@}

clean:	make_kernel_clean
	rm -f .Depend ${OBJS} ${TARGET}
	rm -rf ${OBJDIR}


#######################################################################
#                            Depend processing                        #
#######################################################################

make_kernel: $(KERNEL_DIRS)
	@for d in $(KERNEL_DIRS); \
	do \
	(cd $$d && \
	echo "" && \
	echo "*******************************************************" && \
	echo "* making in $(CURRENT_DIR)/$$d" && \
	echo "*******************************************************" && \
	$(MAKE) $(MFLAGS) "CURRENT_DIR=$(CURRENT_DIR)/$$d" \
	"DEBUG=$(DEBUG)" all); \
	if [ $$? -ne 0 ]; then \
		echo "!!!!!! ERROR !!!!!!" ; exit 1; \
	fi \
	done

make_kernel_clean: $(KERNEL_DIRS)
	@for d in $(KERNEL_DIRS); \
	do \
	(cd $$d && \
	echo "" && \
	echo "*******************************************************" && \
	echo "* making in $(CURRENT_DIR)/$$d" && \
	echo "*******************************************************" && \
	$(MAKE) $(MFLAGS) "CURRENT_DIR=$(CURRENT_DIR)/$$d" \
	"DEBUG=$(DEBUG)" clean); \
	done

#######################################################################
#                            Depend processing                        #
#######################################################################
# .Depend: # ${SRCS}
# 	@if [ ! -d ${OBJDIR} ]; then \
# 		echo ";; mkdir ${OBJDIR}"; mkdir ${OBJDIR}; \
# 	fi
# 	${CC} ${CCOPT} ${INC_DIR} -MM ${CSRC} ${ASRC} > .Depend
#
# -include .Depend
