#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(PDB2IDL, 1.0, bd512@york.ac.uk)

AC_ARG_VAR(CFLAGS,  [Compile flags])
AC_ARG_VAR(LDFLAGS, [Linker flags])
AC_ARG_VAR(IDLPATH, [Path to IDL external include (should contain idl_export.h)])
AC_ARG_VAR(PACT,    [Path to the Portable Application Toolkit library])

#####################################################################
# Compiler / Linker flags
#####################################################################

if test "$CFLAGS" = ""
then
	CFLAGS="-Wall -fPIC -g"
fi

if test	"$LDFLAGS" = ""
then
	LDFLAGS="-shared"
fi

#####################################################################
# IDL
#####################################################################

if test "$IDLPATH" != ""
then
	# Path set. Check it's correct
	AC_CHECK_FILE($IDLPATH/idl_export.h, ,
	AC_CHECK_FILE($IDLPATH/include/idl_export.h, IDLPATH=$IDLPATH/include, 
	AC_CHECK_FILE($IDLPATH/external/include/idl_export.h, IDLPATH=$IDLPATH/external/include,
	IDLPATH=""))) # If not correct, clear path
else
	echo "IDLPATH variable not set. Searching for library"
fi

## Use location of binary
if test "$IDLPATH" = ""
then
	BINPATH=`which idl`
	if test "$BINPATH" != ""
	then
		# This could be a symlink, so try using ls to expand
		# NOTE: ls -L doesn't always seem to work so using -l and sed
		FULLBINPATH=`ls -l $BINPATH | sed "s/ /\n/g"  | tail -n 1`
		if test "$FULLBINPATH" != ""
		then
			FULLBINPATH=`dirname $FULLBINPATH`
			AC_CHECK_FILE($FULLBINPATH/../external/include/idl_export.h, IDLPATH=$FULLBINPATH/../external/include, )

			if test "$IDLPATH" = ""
			then
				# Try using find
				FINDFILE=`find $FULLBINPATH/../ -name idl_export.h | head -n 1`
				if test "$FINDFILE" != ""
				then
				   FINDFILE=`dirname $FINDFILE`
				   AC_CHECK_FILE($FINDFILE/export.h, IDLPATH=$FINDFILE, )
				fi
			fi
		fi
		
		# try find on BINPATH
		if test "$IDLPATH" = ""
		then
			BINPATH=`dirname $BINPATH`
			FINDFILE=`find $BINPATH/../ -name idl_export.h | head -n 1`
			if test "$FINDFILE" != ""
			then
			   FINDFILE=`dirname $FINDFILE`
			   AC_CHECK_FILE($FINDFILE/idl_export.h, IDLPATH=$FINDFILE, )
			fi
		fi
	else
		echo "checking idl command... no"
	fi
fi

## Try to find IDL using locate
if test "$IDLPATH" = ""
then
	IDLFILE=`locate idl_export.h | head -n 1`
	if test "$IDLFILE" != ""
	then
		IDLFILE=`dirname $IDLFILE`
		AC_CHECK_FILE($IDLFILE/idl_export.h, IDLPATH=$IDLFILE, )
	else
		echo "checking locate... no"
	fi
fi

## Try some known paths
if test "$IDLPATH" = ""
then
	AC_CHECK_FILE(/usr/local/rsi/idl/external/include/idl_export.h, IDLPATH=/usr/local/rsi/idl/external/include,
	AC_CHECK_FILE(/afs/localcell/usr/rsi/idl/external/include/idl_export.h, IDLPATH=/afs/localcell/usr/rsi/idl/external/include,
	AC_CHECK_FILE(/usr/common/graphics/idl/idl/external/include/idl_export.h, IDLPATH=/usr/common/graphics/idl/idl/external/include,
	)))
fi

# Final check to see if path was found
if test "$IDLPATH" = ""
then
	echo "Error - IDL external library not found"	
	exit
fi

#####################################################################
# PACT library
#####################################################################

if test "$PACT" = ""
then
	echo "PACT variable not set. Searching for library..." 
	AC_CHECK_FILE(/usr/include/pdb.h, PACT=/usr/,
	AC_CHECK_FILE(/usr/local/include/pdb.h, PACT=/usr/local/,
	AC_CHECK_FILE($HOME/local/include/pdb.h, PACT=$HOME/local/,
	echo "Error - can't find PACT distribution"
	exit)))
fi

# Check that PACT libraries exist
AC_CHECK_FILES( $PACT/include/pdb.h $PACT/lib/libpdb.a $PACT/lib/libpml.a $PACT/lib/libscore.a ,, [
echo "Error - PACT library required"
exit] )


# Checks for programs.
AC_CHECK_PROGS(CC, gcc cc cl)
AC_SUBST(CC)
LD=$CC
AC_SUBST(LD)

# Only make.config is altered by configure
AC_OUTPUT(make.config)

