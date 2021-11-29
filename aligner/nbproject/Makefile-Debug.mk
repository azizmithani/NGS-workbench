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
	${OBJECTDIR}/aligner.o \
	${OBJECTDIR}/bwa.o \
	${OBJECTDIR}/mapper.o \
	${OBJECTDIR}/sam.o


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
LDLIBSOPTIONS=-L../utility/dist/Debug/GNU-MacOSX -lutility

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libaligner.${CND_DLIB_EXT}
	${CP} ../utility/dist/Debug/GNU-MacOSX/libutility.dylib ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	-install_name_tool -change libutility.dylib @executable_path/libutility.dylib ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libaligner.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libaligner.${CND_DLIB_EXT}: ../utility/dist/Debug/GNU-MacOSX/libutility.dylib

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libaligner.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libaligner.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -dynamiclib -install_name libaligner.${CND_DLIB_EXT} -fPIC

${OBJECTDIR}/aligner.o: aligner.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../utility -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/aligner.o aligner.cpp

${OBJECTDIR}/bwa.o: bwa.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../utility -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/bwa.o bwa.cpp

${OBJECTDIR}/mapper.o: mapper.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../utility -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/mapper.o mapper.cpp

${OBJECTDIR}/sam.o: sam.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../utility -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sam.o sam.cpp

# Subprojects
.build-subprojects:
	cd ../utility && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} -r ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libutility.dylib
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libaligner.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:
	cd ../utility && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
