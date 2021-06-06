import numpy as np
from scipy import interpolate
from spectral_difference import local_max
from common import wrap

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
    """
    For use with the amplitudes of the partials at each time step. Any partial
    amplitude that is a local maximum that is more that K-times its neighbours
    is made to be the average of the neighbours. The neighbours are the
    instantaneous amplitudes of the partials just lower or higher in frequency
    at the same time step.
    A is a matrix of shape (n_partials,N) where N is the number of time steps
    and n_partials the number of paritals. Therefore the local maxima are
    searched for in the columns of A.
    """
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

def interpolate_angular_velocity(Th,B=0,rate=1.):
    """
    Interpolates the angular velocities using linear interpolation to obtain an
    interpolated phase sequence.
    Th is a matrix of size (n_partials,N), where N is the number of time steps
    and n_partials the number of partials. Th should represent phases at every
    sample so that there is no ambiguity in the amount the phase has changed
    between consecutive time-points.
    B is the number of columns at the beginning of Th that should not be
    interpolated. The remaining columns (len(Th)-B) are interpolated with a
    lowered rate so that the final number of columns is round(N/rate).
    Returns a matrix of size (n_partials,ceil(N/rate))
    """
    if B != 0:
        raise NotImplementedError
    n_partials=Th.shape[0]
    N=Th.shape[1]
    deltas=np.zeros((n_partials,N))
    deltas[:,1:]=wrap(np.diff(Th))
    n_orig=np.arange(N)
    n_interp=np.linspace(0,N-1,int(np.round(N/rate)))
    delta_interp=interpolate.interp1d(n_orig,deltas)(n_interp)
    Th_interp=Th[:,0][:,None]+np.cumsum(delta_interp,axis=1)
    return Th_interp

def interpolate_amplitudes(A,B=0,rate=1.):
    """
    Works like interpolate_angular_velocity but for amplitudes. B works the same as above.
    Returns a matrix of size (n_partials,ceil(N/rate))
    """
    if B != 0:
        raise NotImplementedError
    n_partials=A.shape[0]
    N=A.shape[1]
    n_orig=np.arange(N)
    n_interp=np.linspace(0,N-1,int(np.round(N/rate)))
    A_interp=interpolate.interp1d(n_orig,A)(n_interp)
    return A_interp

