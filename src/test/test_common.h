#ifndef COMMON_H
#define COMMON_H 

#include <stdlib.h>
#include <stdio.h>

#define PRINT_ANY(x)\
    printf(_Generic( (x), int : "%d", float : "%f"),x)

#define PRINT_ALL(a_,n) \
    ({ __auto_type a__ = (a_);\
       __auto_type n__ = (n);\
       while (n__--) { PRINT_ANY(*a__++); printf(" "); };\
       printf("\n"); })

#define CHK_EQ(a_,b_,n) \
    ({ __auto_type a__ = (a_);\
       __auto_type b__ = (b_);\
       __auto_type n__ = (n);\
       int eq__ = 1 ;\
       while (n__--) { eq__ &= (*a__++ == *b__++); };\
       eq__; })

static inline long
get_file_length(
    FILE *f)
{
    fseek(f,0,SEEK_END);
    long ret = ftell(f);
    rewind(f);
    return ret;
}

static inline const char *
getenv_default(
    const char *env_name,
    const char *default_value)
{
    const char *ret = getenv(env_name);
    if (!ret) { return default_value; }
    return ret;
}

#endif /* COMMON_H */
