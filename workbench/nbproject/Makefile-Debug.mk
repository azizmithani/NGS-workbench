#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/workbench.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L../aligner/dist/Debug/GNU-MacOSX -laligner -L../utility/dist/Debug/GNU-MacOSX -lutility

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/workbench
	${CP} ../aligner/dist/Debug/GNU-MacOSX/libaligner.dylib ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	-install_name_tool -change libaligner.dylib @executable_path/libaligner.dylib ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/workbench
	${CP} ../utility/dist/Debug/GNU-MacOSX/libutility.dylib ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	-install_name_tool -change libutility.dylib @executable_path/libutility.dylib ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/workbench

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/workbench: ../aligner/dist/Debug/GNU-MacOSX/libaligner.dylib

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/workbench: ../utility/dist/Debug/GNU-MacOSX/libutility.dylib

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/workbench: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/workbench ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/workbench.o: workbench.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../aligner -I../utility -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/workbench.o workbench.cpp

# Subprojects
.build-subprojects:
	cd ../aligner && ${MAKE}  -f Makefile CONF=Debug
	cd ../utility && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} -r ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libaligner.dylib ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libutility.dylib
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/workbench

# Subprojects
.clean-subprojects:
	cd ../aligner && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../utility && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
