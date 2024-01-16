# Time-Frequency Representations

def reassigned_stft(x,n,N,h,dh):
    """
    Compute the reassigned STFT of x.

    x in R^len(x) or C^len(x): the signal to compute on
    n in Z^len(n): the sample times
    N a positive integer: the size of DFT transform. This determines the
                          frequency resolution of the output. The best
                          performance is obtained when this is a power of 2.
    h in R^len(x): the window applied to each N-sized subsequence of x. Only
                   windows of length less than or equal to N are supported.
    dh in R^len(x): the derivative w.r.t to n_w, the sample times of the window.
                    The sample times of the window are always centred around each n,
                    i.e., [len(h)/2,...,len(h)/2]

    Returns (n_r,v_r,X) with
    n_r in R^len(n)xN the reassigned sample times
    v_r in R^len(n)xN the reassigned normalized frequencies
    X   in C^len(n)XN the spectrogram
    """

    
