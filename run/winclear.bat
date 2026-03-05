rmdir /s /q .vs
rmdir /s /q CMakeFiles
rmdir /s /q src
rmdir /s /q x64

rmdir /s /q Debug
rmdir /s /q Release
rmdir /s /q MinSizeRel
rmdir /s /q RelWithDebInfo


del ALL_BUILD.*
del ZERO_CHECK.*
del cmake_install.*
del CMakeCache.txt
del VM.sln