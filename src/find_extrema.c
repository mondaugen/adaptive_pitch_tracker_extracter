
/*
def discount_local_max(x,rate,min_thresh=0):
    lmaxs=local_max(x)
    thresh=np.zeros_like(x)
    filtered_maxs=[]
    for n_max in lmaxs:
        if (x[n_max] > min_thresh) and (x[n_max] > thresh[n_max]):
            s=np.zeros_like(x)
            s[n_max]=x[n_max]
            thresh+=signal.lfilter([1],[1,-rate],s)
            filtered_maxs.append(n_max)
    return (np.array(filtered_maxs),thresh)
*/

/* if 'right', then a point >= to a point to its right is a candidate for a
local maximum */
#define sided_max_right(x_1,x0,x1) ((x0>x_1)&&(x0>=x1))

/* if 'left', then a point >= to a point to its left is a candidate for a local
maximum */
#define sided_max_left(x_1,x0,x1) ((x0>=x_1)&&(x0>x1))

/* if 'none', then a point must be > than both points */
#define sided_max_none(x_1,x0,x1) ((x0>x_1)&&(x0>x1))

/* if 'both', then a point only needs to be >= than both points */
#define sided_max_both(x_1,x0,x1) ((x0>=x_1)&&(x0>=x1))

#define sided_max(a,b,c,t)\
    ({\
    int r = 0;\
    switch (t) {\
    case local_max_type_right:\
        r = sided_max_right(a,b,c);\
        break;\
    case local_max_type_left:\
        r = sided_max_left(a,b,c);\
        break;\
    case local_max_type_none:\
        r = sided_max_none(a,b,c);\
        break;\
    case local_max_type_both:\
        r = sided_max_both(a,b,c);\
        break;\
    }\
    r })

unsigned int *
local_max_f32(const float *x, /* find local maxima in this array */
              unsigned int length, /* length of array to search */
              unsigned int *n_max, /* number of local maxima found */ 
              enum local_max_type type /* maxima search type */)
{
    unsigned int n, cnt = 0;
    /* First count the number of local maxima */
    for (n = 1; n < (length-1); n++) {
        if (sided_max(x[n-1],x[n],x[n+1],type)) { cnt += 1; }
    } 
    /* Then allocate the array to hold them */
    unsigned int *ret = malloc(sizeof(unsigned int [cnt]));
    if (!ret) { return NULL; }
    *n_max = cnt;
    cnt = 0;
    /* Now store the maxima */
    for (n = 1; n < (length-1); n++) {
        if (sided_max(x[n-1],x[n],x[n+1],type)) {
            ret[cnt] = n;
            cnt += 1;
        }
    } 
    return ret;
}

/*
When a local maximum is encountered, it is compared with a threshold
function (which is initially 0). If it is greater than the function at its
point then the threshold function has this max convolved with an exponential
decay summed into it, the maximum is recorded, and the algorithm proceeds.
If it is not then this maximum is discarded and the algorithm proceeds.
This returns the filtered maxima and the threshold function
The value must be over min_thresh to be accepted.
*/
unsigned int *
discount_local_max_f32(
    /* find local maxima in this array */
    const float *x,
    /* length of array to search */
    unsigned int length, 
    /* number of local maxima found */ 
    unsigned int *n_max, 
    /* maxima search type */
    enum local_max_type type, 
    /* discount rate */
    float rate, 
    /* threshold above which maximum is accepted */
    float min_thresh, 
    /* if not NULL, must be an array of length "length" and will
    contain the threshold function after the function returns */
    float *treshold)
{
    unsigned int *lmaxs = NULL, lmaxs_len,
                 *filtered_maxs = NULL, filt_maxs_len = 0,
                 n_x, n_maxs = 0;
    float cur_thresh = 0, th_x0;
    lmaxs = local_max_f32(x,length,&lmaxs_len,type);
    if (!lmaxs) { goto fail; }
    /* Count the number of maxima and store the threshold function (if not NULL) */
    for (n_x = 0; n_x < length; n_x++) {
        /* the x[n] coefficient of the threshold updating filter */
        th_x0 = 0;
        if ((n_x == lmaxs[n_maxs])
            &&(x[n_x]>min_thresh)
            &&(x[n_x]>(cur_thresh*rate))) {
            filt_maxs_len += 1;
            th_x0 = x[n_x];
            n_maxs++;
        }
        cur_thresh = th_x0 + rate * cur_thresh;
        if (thresh) { thresh[n_x] = cur_thresh; }
    }
    /* allocate the memory necessary for the filtered maxima */
    filtered_maxs = malloc(sizeof(unsigned int [n_maxs]));
    if (!filtered_maxs) { goto fail; }
    cur_thresh = 0;
    n_maxs = 0;
    for (n_x = 0; n_x < length; n_x++) {
        /* the x[n] coefficient of the threshold updating filter */
        th_x0 = 0;
        if ((n_x == lmaxs[n_maxs])
            &&(x[n_x]>min_thresh)
            &&(x[n_x]>(cur_thresh*rate))) {
            filt_maxs_len += 1;
            th_x0 = x[n_x];
            n_maxs++;
        }
        cur_thresh = th_x0 + rate * cur_thresh;
    }
fail:
    if (lmaxs) { free(lmaxs); }
    return filtered_maxs;
}

                 
    lmaxs=local_max(x)
    thresh=np.zeros_like(x)
    filtered_maxs=[]
    for n_max in lmaxs:
        if (x[n_max] > min_thresh) and (x[n_max] > thresh[n_max]):
            s=np.zeros_like(x)
            s[n_max]=x[n_max]
            thresh+=signal.lfilter([1],[1,-rate],s)
            filtered_maxs.append(n_max)
    return (np.array(filtered_maxs),thresh)
}
    
