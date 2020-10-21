struct adder {
    float v;
};

struct adder_proc {
    unsigned int N;
    float *a;
    float *b;
};

void adder_add(struct adder_proc *p) {
    unsigned int N = p->N;
    while (N--) {
        p->a[N] += p->b[N];
    }
}
