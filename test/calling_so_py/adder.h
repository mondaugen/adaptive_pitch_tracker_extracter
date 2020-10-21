#ifndef ADDER_H
#define ADDER_H

struct adder;

struct adder_proc {
    unsigned int N;
    float *a;
    float *b;
};

void adder_add(struct adder_proc *p);

#endif /* ADDER_H */
