import numpy as np
from polyeval import polyeval

RECIP_2PI=1./(2.*np.pi)
TWO_PI=2.*np.pi

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
        self.theta_k0[:]=theta_k0
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
                                + 0.5*self.H(omega_k1 - omega_k0)))
        # determine coefficients of smoothest cubic polynomial
        alpha_beta_M=self.alpha_beta @ np.vstack(
            theta_k1 - self.theta_k0 - self.H * self.omega_k0 + TWO_PI*M,
            omega_k1 - self.omega_k0
        )
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
        


