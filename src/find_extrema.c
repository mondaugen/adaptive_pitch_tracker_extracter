
/*
def discount_local_max(x,rate,min_thresh=0):
    """
    When a local maximum is encountered, it is compared with a threshold
    function (which is initially 0). If it is greater than the function at its
    point then the threshold function has this max convolved with an exponential
    decay summed into it, the maximum is recorded, and the algorithm proceeds.
    If it is not then this maximum is discarded and the algorithm proceeds.
    This returns the filtered maxima and the threshold function
    The value must be over min_thresh to be accepted.
    """
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

