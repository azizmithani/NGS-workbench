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
	${OBJECTDIR}/base_distribution.o \
	${OBJECTDIR}/dna.o \
	${OBJECTDIR}/gff.o \
	${OBJECTDIR}/gff_entry.o \
	${OBJECTDIR}/position_base_distribution.o \
	${OBJECTDIR}/region.o \
	${OBJECTDIR}/snp.o \
	${OBJECTDIR}/utility.o


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
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libutility.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libutility.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libutility.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -dynamiclib -install_name libutility.${CND_DLIB_EXT} -fPIC

${OBJECTDIR}/base_distribution.o: base_distribution.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/base_distribution.o base_distribution.cpp

${OBJECTDIR}/dna.o: dna.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/dna.o dna.cpp

${OBJECTDIR}/gff.o: gff.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/gff.o gff.cpp

${OBJECTDIR}/gff_entry.o: gff_entry.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/gff_entry.o gff_entry.cpp

${OBJECTDIR}/position_base_distribution.o: position_base_distribution.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/position_base_distribution.o position_base_distribution.cpp

${OBJECTDIR}/region.o: region.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/region.o region.cpp

${OBJECTDIR}/snp.o: snp.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/snp.o snp.cpp

${OBJECTDIR}/utility.o: utility.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/utility.o utility.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
