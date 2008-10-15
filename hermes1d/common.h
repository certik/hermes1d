#ifndef common_h
#define common_h

#include <stdlib.h>

inline void error(const char *str)
{
    printf("ERROR: %s\n", str);
    exit(1);
}

#endif
