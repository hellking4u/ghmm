# Microsoft Developer Studio Project File - Name="ghmmpp" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=GHMMPP - WIN32 RELEASE
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ghmmpp.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ghmmpp.mak" CFG="GHMMPP - WIN32 RELEASE"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ghmmpp - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe
# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "..\..\\" /I "..\..\..\gsl-1.0" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x407 /d "NDEBUG"
# ADD RSC /l 0x407 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
# Begin Target

# Name "ghmmpp - Win32 Release"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_AbstractModel.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_ContinuousModel.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_DiscreteModel.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Document.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_DoubleMatrix.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_DoubleVector.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Emission.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_IntVector.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Sequence.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Sequences.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_State.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Toolkit.cpp"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Transition.cpp"
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE="..\..\ghmm++\begin_code.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\close_code.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_AbstractModel.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_ContinuousModel.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_DiscreteModel.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Document.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_DoubleMatrix.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_DoubleVector.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Emission.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_IntVector.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Sequence.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Sequences.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_State.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Transition.h"
# End Source File
# Begin Source File

SOURCE="..\..\ghmm++\GHMM_Types.h"
# End Source File
# End Group
# End Target
# End Project
