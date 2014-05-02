cmake ^
-D CMAKE_TOOLCHAIN_FILE=%PRODDL_SRC%\config\win\toolchain.cmake ^
-G "Visual Studio 10" ^
..\deps_src\rmf

REM "C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\devenv" PRODDL.sln /build Release /OUT devenv.build.log /Project ALL_BUILD.vcxproj