#----------------------------------------------------------#
# Upper level makefile, determines OS and calls makefile   #
# for this specific OS                                     #
#----------------------------------------------------------#

CPU = $(shell uname -s)


ifeq ($(CPU), AIX)
all:
	gmake -f efitxs_aix.mak
clean:
	gmake -f efitxs_aix.mak clean
endif


ifeq ($(CPU), Linux)
all:
	gmake -f efitxs_linux.mak
clean:
	gmake -f efitxs_linux.mak clean
endif
