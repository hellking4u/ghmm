# Microsoft Developer Studio Project File - Name="ghmm" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=ghmm - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ghmm.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ghmm.mak" CFG="ghmm - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ghmm - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "ghmm - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "ghmm - Win32 Release"

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

!ELSEIF  "$(CFG)" == "ghmm - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "..\..\\" /I "..\..\..\gsl-1.0" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x407 /d "_DEBUG"
# ADD RSC /l 0x407 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "ghmm - Win32 Release"
# Name "ghmm - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\ghmm\cluster.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\foba.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\gauss_tail.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\matrix.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\mes.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\model.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\mprintf.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\randvar.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\reestimate.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\rng.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\root_finder.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\scanner.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\scluster.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sequence.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sfoba.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sgenerate.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\smap_classify.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\smixturehmm.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\smodel.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sreestimate.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sviterbi.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\vector.c
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\viterbi.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\ghmm\cluster.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\const.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\foba.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\gauss_tail.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\matrix.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\mes.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\model.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\mprintf.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\randvar.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\reestimate.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\rng.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\root_finder.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\scanner.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\scluster.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sequence.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sfoba.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sgenerate.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\smap_classify.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\smixturehmm.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\smodel.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sreestimate.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\sviterbi.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\vector.h
# End Source File
# Begin Source File

SOURCE=..\..\ghmm\viterbi.h
# End Source File
# End Group
# End Target
# End Project
