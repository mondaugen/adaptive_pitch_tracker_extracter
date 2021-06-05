import numpy as np
from spectral_difference import local_max

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
    
def column_max(X,one_sided_max='right',K=1,return_bool=False):
    X_=X.T.flatten()
    max_mask=local_max(X_,one_sided_max=one_sided_max,K=K,return_bool=True)
    # Maxima cannot be at the beginning or end of columns
    max_mask[0::X.shape[0]]=False
    max_mask[X.shape[0]-1::X.shape[0]]=False
    max_mask_mat=max_mask.reshape(X.shape[::-1]).T
    if return_bool:
        return max_mask_mat
    return np.where(max_mask_mat)

def squish_anomalies(A,K=2):
    max_mask=column_max(A,one_sided_max='none',K=K,return_bool=True)
    fill_vals=np.zeros_like(A)
    # Make shifted masks to get the values one above and one below the local max
    # values
    mask_shifted_up=np.zeros_like(max_mask)
    mask_shifted_up[:-1]=max_mask[1:]
    mask_shifted_down=np.zeros_like(max_mask)
    mask_shifted_down[1:]=max_mask[:-1]
    # Find the average of values one above and one below the local max values
    fill_vals[max_mask]=A[mask_shifted_up]
    fill_vals[max_mask]+=A[mask_shifted_down]
    fill_vals*=0.5
    A[max_mask]=fill_vals[max_mask]
    return A




