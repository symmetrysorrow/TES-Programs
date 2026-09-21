@ECHO OFF
SET LIBDIR=%ELMER_HOME%/bin
SET INCLUDE=%ELMER_HOME%/share/elmersolver/include

SET FC="%ELMER_HOME%/stripped_gfortran/bin/"

SET cmd=%FC% %*  -fallow-argument-mismatch  -DELMER_BROKEN_MPI_IN_PLACE -DCONTIG= -DMINGW32 -DWIN32 -DHAVE_EXECUTECOMMANDLINE -DUSE_ISO_C_BINDINGS -DUSE_ARPACK   -shared -I"%INCLUDE%" -L"%LIBDIR%" -shared -lelmersolver
echo %cmd%
%cmd%
