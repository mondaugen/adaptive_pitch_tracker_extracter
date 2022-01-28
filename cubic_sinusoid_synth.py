import numpy as np
from polyeval import polyeval
from scipy import interpolate
from common import wrap

RECIP_2PI=1./(2.*np.pi)
TWO_PI=2.*np.pi
j=complex('j')

def cubic_phase_poly_interp(H,theta,omega):
    """
    Computes the polynomial coefficients for the smoothest cubic phase
    polynomial splines and evaluates these polynomials. Unlike
    cubic_sinusoid_synth.process, this can process theta and omega as matrices,
    and therefore all the splines in a batch.
    H (an integer) is the number of samples between phase measurements.
    theta is a matrix of size (F,n_partials) representing the phases of the
    n_partials sinusoids at each breakpoint.
    omega is a matrix of size (F,n_partials) representing the angular velocities
    at each breakpoint.
    The result is a matrix of size (n_partials,(F-1)*H+1) containing the phases
    of the n_partials sinusoids evaluated at sample 0, 1, ...., (F-1)*H+1.
    """
    F=theta.shape[0]
    n_partials=theta.shape[1]
    assert (F == omega.shape[0]) and (n_partials == omega.shape[1])
    # matrix that is used to determine phase polynomial coefficients
    alpha_beta=np.array([
        [-2./(H*H*H), 1./(H*H)],
        [3./(H*H), -1./H]
    ])
    theta_k0=theta[:-1,:]
    theta_k1=theta[1:,:]
    omega_k0=omega[:-1,:]
    omega_k1=omega[1:,:]
    # value used to determine smoothest interpolating phase polynomial
    M=np.round(RECIP_2PI*(theta_k0+omega_k0*H - theta_k1 
                            + 0.5*H*(omega_k1 - omega_k0)))
    # coefficients of smoothest cubic polynomial
    omega_theta_M=np.stack((
        theta_k1 - theta_k0 - H * omega_k0 + TWO_PI*M,
        omega_k1 - omega_k0
    )).transpose(1,0,2)
    pcoefs=np.zeros((F-1,4,n_partials))
    pcoefs[:,:2,:]=alpha_beta @ omega_theta_M
    pcoefs[:,2,:]=omega_k0
    pcoefs[:,3,:]=theta_k0
    n=np.arange(H)
    Th=np.zeros((n_partials,(F-1)*H+1))
    Th[:,:-1]=np.hstack(polyeval(np.hstack(pcoefs),n).reshape((F-1,n_partials,H)))
    Th[:,-1]=theta[-1]
    return Th

def single_frame_quadratic_phase_poly_interp(th,o0,o1,H):
    """
    Compute a single frame of polynomial interpolation
    theta has size n_partials, the initial phase (time h=0)
    o0 has size n_partials, the initial frequencies (time h=0)
    o1 has size n_partials, the frequencies at time h=H
    Returns a matrix of size(n_partials,H+1)
    """
    a=(o1-o0)/(2.*H)
    b=o0
    c=th
    pcoeffs=np.vstack((a,b,c))
    h=np.arange(H+1)
    phs=polyeval(pcoeffs,h)
    return phs

def quadratic_phase_poly_interp(H,theta,omega):
    """
    Computes the polynomial coefficients for the piece-wise quadratic splines
    and evaluates these polynomials.
    H (an integer) is the number of samples between phase measurements.
    theta is a vector of size n_partials containing the initial phases of the
    n_partials sinusoids.
    omega is a matrix of size (F,n_partials) representing the angular velocities
    at each breakpoint
    The result is a matrix of size (n_partials,(F-1)*H+1) containing the phases
    of the n_partials sinusoids evaluated at sample 0, 1, ...., (F-1)*H+1.
    """
    F=omega.shape[0]
    n_partials=theta.shape[0]
    assert n_partials == omega.shape[1]
    cur_theta=theta
    ret=np.zeros((n_partials,(F-1)*H+1))
    for f in np.arange(F-1):
        start_sample=H*f
        end_sample=H*(f+1)+1
        phs=single_frame_quadratic_phase_poly_interp(cur_theta,omega[f,:],omega[f+1,:],H)
        phs=wrap(phs)
        ret[:,start_sample:end_sample]=phs
        cur_theta=ret[:,end_sample-1]
    return ret


def linear_amplitude_interp(H,a):
    """
    Takes a matrix "a" of shape (F,n_partials) containing F amplitude estimates
    (taken every H samples) for n_partials partial tracks.
    The result is a matrix of size (n_partials,(F-1)*H+1) containing the
    amplitudes of the n_partials for every sample by interpolating linearly
    between the amplitude estimates.
    """
    F=a.shape[0]
    n=np.arange(F)*H
    ni=np.linspace(0,n[-1],(F-1)*H+1)
    A=interpolate.interp1d(n,a,axis=0)(ni)
    return A.T
    

class cubic_sinusoid_synth:
    def __init__(self,n_partials,H):
        """
        Create a new cubic_sinusoid_synth.
        n_partials is the maximum number of partial tracks that can be
        synthesized in parallel (it will store this many phases between
        synthesis calls).
        H is the hop size in samples.
        """
        assert(H != 0)
        self.H=H
        # last phase (in radians)
        self.theta_k0=np.zeros(n_partials)
        # last angular velocity (in radians/sample)
        self.omega_k0=np.zeros(n_partials)
        # last amplitude (in linear amplitude)
        self.a_k0=np.zeros(n_partials)
        # matrix that is used to determine phase polynomial coefficients
        self.alpha_beta=np.array([
            [3./(H*H), -1./H],
            [-2./(H*H*H), 1./(H*H)]
        ])
        # the sample numbers that are used for polynomial evaluation
        self.n=np.arange(H)
        # The reciprocal of H, used for amplitude interpolation.
        self.recip_H=1./self.H
    def set_theta_k0(self,theta_k0):
        self.theta_k0[:]=np.fmod(theta_k0,2*np.pi)
    def set_omega_k0(self,omega_k0):
        self.omega_k0[:]=omega_k0
    def set_a_k0(self,a_k0):
        self.a_k0[:]=a_k0
    def process(self,theta_k1,omega_k1,a_k1):
        """
        Compute H samples of n_partials tracks.
        theta_k1 is the phase after H samples.
        omega_k1 is the angular velocity after H samples.
        a_k1 is the amplitude after H samples.

        The synthesis is done as in McAulay and Quatieri (1986). The phase is
        interpolated to maximize smoothness while obeying the phases and angular
        velocities at the end points, and the amplitude is interpolated
        linearly.

        After the synthesis is finished, the {theta,omega,a}_k1 are stored in
        {theta,omega,a}_k0 to serve as the starting values for the next call to
        process.

        Returns: Matricies Th, and A. 
        Th is a matrix of size (n_partials,H) containing the phases of the
        n_partials sinusoids evaluated at sample 0, 1, ... H-1. These can then
        be used in a lookup call as exp(j*Th) or cos(Th).
        A is a matrix of size (n_partials, H) containing the amplitude values
        evaluated at sample 0, 1, ... H-1. These can be used as-is or in a
        lookup call as exp(A) (if a_k{0,1} are log amplitude values, for
        example).
        """

        ## Phase computation

        # value used to determine smoothest interpolating phase polynomial
        M=np.round(RECIP_2PI*(self.theta_k0+self.omega_k0*self.H - theta_k1 
                                + 0.5*self.H*(omega_k1 - self.omega_k0)))
        # determine coefficients of smoothest cubic polynomial
        alpha_beta_M=self.alpha_beta @ np.vstack((
            theta_k1 - self.theta_k0 - self.H * self.omega_k0 + TWO_PI*M,
            omega_k1 - self.omega_k0
        ))
        # evaluate the polynomial
        pcoefs = np.vstack((np.flipud(alpha_beta_M), self.omega_k0,
        self.theta_k0))
        Th=polyeval(pcoefs,self.n)

        ## Amplitude computation

        a_diff = a_k1 - self.a_k0
        A = self.a_k0[:,None] + np.multiply.outer(a_diff*self.recip_H,self.n)

        # update parameters
        self.set_theta_k0(theta_k1)
        self.set_omega_k0(omega_k1)
        self.set_a_k0(a_k1)

        return (Th,A)
        

def process_partial_tracks(H,theta,omega,a):
    """
    Computes the partial tracks using the theta, omega and a as breakpoints.
    H is the number of samples between breakpoints.
    theta is a matrix of size (F,n_partials) representing the phases of the
    n_partials sinusoids at each breakpoint.
    omega is a matrix of size (F,n_partials) representing the angular velocities
    at each breakpoint.
    a is a matrix of size (F,n_partials) representing the amplitudes at each
    breakpoint.
    Note that theta, omega and a must be the same size.
    Returns: matricies Th, and A.
    Th is a matrix of size (n_partials,H*(F-1)), representing the phases of the n_partials.
    A is matrix of size (n_partials,H*(F-1)), containing the amplitudes of the n_partials.
    """
    assert(theta.shape == omega.shape == a.shape)
    n_partials=theta.shape[1]
    F=theta.shape[0]
    part_synth=cubic_sinusoid_synth(n_partials,H)
    part_synth.set_theta_k0(theta[0])
    part_synth.set_omega_k0(omega[0])
    part_synth.set_a_k0(a[0])
    Th=cubic_phase_poly_interp(H,theta,omega)[:,:H*(F-1)]
    A=linear_amplitude_interp(H,a)[:,:H*(F-1)]
    return (Th,A)

def synth_partial_tracks(Th,A,th_mode='cos',a_mode='linear',combine=True,
                         phase_units='radians'):
    """
    Takes matrix Th of size (n_partials,N) of phases and matrix A of size
    (n_partials,N) of amplitudes and uses these to synthesize the partial
    tracks.
    th_mode can be 'cos' or 'cexp' or a function taking a matrix as argument and
    returning a matrix of the same size. 'cos' uses the phases as the argument
    of the cosine function, 'cexp' uses the phases as the argument of the
    complex exponential function.
    a_mode can be 'linear', 'dB' or a function taking a matrix as argument and
    returning a matrix of the same size. Linear just uses the amplitude values
    directly, dB converts from decibels to amplitude, using the function
    10^(A/20).
    combine, if true multiplies each sinusoid track by the corresponding
    amplitude track, sums these down the columns and returns an array of length
    N, otherwise it returns matrices (Th_synth,A_synth) which represent the
    sinusoid and amplitude tracks.
    phase_units can be 'radians' or 'normalized'. If 'radians' the phase values
    are used as-is to the arguments of the functions, otherwise the cos(x)
    function is instead cos(2*pi*x) and exp(j*x) is exp(j*2*pi*x). The argument
    to a function passed into th_mode is not affected: the user must provide a
    function accepting normalized phase values in this case.
    """
    # Determine phase processing function
    phase_mul=2*np.pi if phase_units == 'normalized' else 1.
    th_fun=None
    if th_mode == 'cos':
        th_fun=lambda x: np.cos(phase_mul*x)
    elif th_mode == 'cexp':
        th_fun=lambda x: np.exp(j*phase_mul*x)
    else:
        th_fun=th_mode
    # Determine amplitude processing function
    a_fun=None
    if a_mode == 'linear':
        a_fun=lambda x: x
    elif a_mode == 'dB':
        a_fun=lambda x: np.power(10.,x/20)
    # Synthesize tracks
    Th_synth=th_fun(Th)
    A_synth=a_fun(A)
    if combine == False:
        return (Th_synth,A_synth)
    Th_synth*=A_synth
    Th_synth=np.sum(Th_synth,axis=0)
    return Th_synth
