#ifndef COMMON_H
#define COMMON_H 

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

#endif /* COMMON_H */
