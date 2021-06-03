import numpy as np

def common_partial_shape(A,weighted='equal'):
    """
    Takes matrix A of shape (n_partials,N) representing the amplitudes of
    n_partials partials at N consecutive time points. Each partial is normalized
    and the average 'shape' is found by averaging the normalized curves.
    Here normalized means the curve's range is 0-1.
    If weighted='power' then the weighted average is found by weighting each
    curve by its relative power.
    """
    maxes=A.max(axis=1)
    mins=A.min(axis=1)
    diffs=maxes-mins
    ret=A-mins[:,None]
    ret/=diffs[:,None]
    if weighted=='equal':
        return np.mean(ret,axis=0)
    elif weighted=='power':
        # A is real so the squared modulus is A*A
        powers=np.sqrt(np.mean(A*A,axis=1))
        powers/=np.sum(powers)
        return np.sum(ret*powers[:,None],axis=0)
    else:
        raise ValueError('Bad value for "weighted": %s' % (weighted,))

def apply_partial_shape(A,shape):
    """
    Takes a partial shape, for example as determined using common_partial_shape,
    and applies it to each partial in A as follows: (see common_partial_shape
    for the shape of A)
        - Find the minimum and maximum of each partial in A
        - Scale shape to have that range for each partial
    """
    mins=A.min(axis=1)
    maxes=A.max(axis=1)
    diffs=maxes-mins
    A_new=np.multiply.outer(diffs,shape)+mins[:,None]
    return A_new
    
