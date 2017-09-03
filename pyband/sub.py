import numpy as np
from . import band

## Nonlinear Solver
def bound_val(cold, deqfun, order, h, tol, itmax):  

    ( N, NJ )=cold.shape

    for iter in range(1,itmax):
        
        ( dcdx , d2cdx2 )=interpolate( cold, h, order)
        (sma, smb, smd, smg)=deqfun(cold, dcdx, d2cdx2)
        (ABD,G)=buildMAT(sma, smb, smd, smg, h)

        delc= band.abdgxy(ABD,G)

        error=np.amax(np.absolute(delc))
        #print "iter, error = %i, %g" % (iter,error)

        cold=cold+delc

        if error < tol:
            return cold

    print 'The program did not converge!!'
    return


## Discretize Derivatives

def interpolate( c, h, order ):

    ( N, NJ )=c.shape
    
    cE=np.concatenate( (         c[:,1:NJ], np.zeros((N,1)) ), axis=1 )
    cW=np.concatenate( (   np.zeros((N,1)),     c[:,0:NJ-1] ), axis=1 )

    dcdx=(cE-cW)/2.0/h;
    dcdx[:, 0  ] = ( - 3.0*c[:, 0  ] + 4.0*c[:,1   ] - c[:,2   ] ) /2.0 /h;
    dcdx[:,NJ-1] = (   3.0*c[:,NJ-1] - 4.0*c[:,NJ-2] + c[:,NJ-3] ) /2.0 /h;

    if ( order >= 2 ):

        d2cdx2= (cE + cW - 2*c)/h**2;
        d2cdx2[:,0] = np.zeros((1,N));
        d2cdx2[:,NJ-1] = np.zeros((1,N));

    return ( dcdx , d2cdx2 )
    

# ### Build Matrix

def buildMAT(sma, smb, smd, smg, h):
    
    ( N, NJ )=sma.shape[0:2]
    
    sma = np.transpose(sma, (0, 2, 1))
    smb = np.transpose(smb, (0, 2, 1))
    smd = np.transpose(smd, (0, 2, 1))
    
    A = sma-h/2.0*smb
    B = -2.0*sma+h**2*smd
    D = sma+h/2.0*smb
    G = h**2*smg
    
    B[:,:,0] = h*smd[:,:,0]-1.5*smb[:,:,0]
    D[:,:,0] = 2.0*smb[:,:,0]
    G[:,0]=h*smg[:,0]
    X = -0.5*smb[:,:,0]
    
    A[:,:,NJ-1]=-2.0*smb[:,:,NJ-1]
    B[:,:,NJ-1]=h*smd[:,:,NJ-1]+1.5*smb[:,:,NJ-1]
    G[:,NJ-1]=h*smg[:,NJ-1]
    Y=0.5*smb[:,:,NJ-1]
    
    ABD = np.concatenate((A, B, D), axis=1)
    BC1 = np.concatenate((B[:,:,0] , D[:,:,0] , X), axis=1)
    BC2 = np.concatenate((Y , A[:,:,NJ-1] , B[:,:,NJ-1]), axis=1)
    ABD[:,:,0] = BC1
    ABD[:,:,NJ-1] = BC2
    
    return ABD, G
