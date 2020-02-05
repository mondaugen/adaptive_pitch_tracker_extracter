#ifndef INDEX_POINTS_H
#define INDEX_POINTS_H 

struct index_points {
    const unsigned int length;
    unsigned int indices[];
};

static inline struct index_points *
index_points_new(unsigned int length);

void
index_points_free(struct index_points *x);

#endif /* INDEX_POINTS_H */
