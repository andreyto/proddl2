cmake ^
-D DEVEL_MODE=1 ^
-D PRODDL_TARGET_ENV=win ^
-D CMAKE_TOOLCHAIN_FILE=%PRODDL_SRC%\config\win\toolchain.cmake ^
-G "Visual Studio 10" ^
%PRODDL_SRC%

cmake --build . --config Release
cmake --build . --target INSTALL.vcxproj --config Release
REM cmake --build . --target PACKAGE.vcxproj --config Release
REM "C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\devenv" PRODDL.sln /build Release /OUT devenv.build.log /Project INSTALL.vcxproj
