import numpy as np
# Gradient adaptive lattice filter as described in Adaptive Signal Processing by
# Alexander p. 106

def ngal(x,P,alpha=0.01,beta=0.01,normalize=True):
    # Produces P refelection coefficients
    N=len(x)
    # allocate memory
    # All have N+1 columns because we need 1 value in the past at each step
    # forward errors, N rows, P+1 columns
    Ef=np.zeros((N+1,P+1))
    # backward errors, N rows, P+1 columns
    Eb=np.zeros((N+1,P+1))
    # adaptivity coefficients
    D=np.zeros((N+1,P))
    D[0,:]=1
    # reflection coefficients
    K=np.zeros((N+1,P))
    Ef[1:N+1,0]=x
    Eb[1:N+1,0]=x
    for n in range(1,N+1):
        for p in range(1,P+1):
            Ef[n,p] = Ef[n,p-1] + K[n-1,p-1] * Eb[n-1,p-1]
            Eb[n,p] = Eb[n-1,p-1] + K[n-1,p-1] * Ef[n,p-1]
            if normalize:
                D[n,p-1] = (1-alpha)*D[n-1,p-1] +  alpha * (Ef[n,p]**2 + Eb[n,p]**2)
            else:
                D[n,p-1] = 1
            K[n,p-1] = K[n-1,p-1] - beta*(Ef[n,p]*Eb[n-1,p-1]
                + Eb[n,p]*Ef[n,p-1]) / D[n,p-1]
    return Ef,Eb,D,K

