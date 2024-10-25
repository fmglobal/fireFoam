#include "stdio.h"

#ifdef CMAKE_VERSION_PRINT
#include "git_version.H"
#else
#include "gitVersion.H"
#endif

if (Pstream::master()){

#ifdef CMAKE_VERSION_PRINT
printf("FireFOAM Build Version:  ");
printf("%10s\n",BUILD_GIT_REVISION);
printf("FireFOAM Build Time Stamp:  ");
printf("%10s\n",BUILD_TIME_STAMP);
#else
printf("FireFOAM Build Version:  ");
printf("%10s\n",build_git_sha);
printf("FireFOAM Build Time Stamp:  ");
printf("%10s\n",build_git_time);
#endif
printf("****************************************************\n");
}
