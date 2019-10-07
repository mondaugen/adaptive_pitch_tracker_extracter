#include <stddef.h>
#include <stdio.h>

struct A {
    const float *f;
    unsigned int n;
    int b;
    float * const g;
};

void
print_A(struct A *a)
{
    printf("f: %p\n"
           "n: %u\n"
           "b: %d\n"
           "g: %p\n",
           a->f,a->n,a->b,a->g);
}

static float g = 123;

struct A
get_pre_init_A(void)
{
    static struct A a = { .g = &g };
    return a;
}

int main (void)
{
    float f[] = {1,3,2,4};
    struct A a = {f,4,6,f},
             b = get_pre_init_A(),
             c = {.g = &g};
    a.f = NULL;
    // Would fail
    // a.g = NULL;
    print_A(&b);
    print_A(&c);
    return 0;
}
