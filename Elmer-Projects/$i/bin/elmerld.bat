@ECHO OFF
SET LIBDIR=%ELMER_HOME%/bin
SET INCLUDE=%ELMER_HOME%/share/elmersolver/include

SET LD="%ELMER_HOME%/stripped_gfortran/bin/"

REM SET CMD=%LD% -shared %* -L"%LIBDIR%" -L"%ELMER_HOME%/bin" -lelmersolver
SET CMD=%LD%  -fallow-argument-mismatch  -shared %* -L"D:/Github/TES-Programs/Elmer-Projects/$i/lib/elmersolver" -L"%LIBDIR%" -lelmersolver
echo %cmd%
%cmd%
