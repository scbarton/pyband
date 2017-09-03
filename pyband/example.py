from numpy import *
from pyband.sub import *

def fillmat(cold,dcdx,d2cdx2):  
 
    #________first column refers to equation
    #________second column refers to position
    #________third column refers to species   
 
    ( N, NJ )=cold.shape
    
    sma=zeros((N,NJ,N))
    smb=zeros((N,NJ,N))
    smd=zeros((N,NJ,N))
    smg=zeros((N,NJ))
 
    sma[0,:,0]=  ones((1,NJ))
    smd[0,:,0]= -cold[1,:]
    smd[0,:,1]= -cold[0,:]
    smg[0,:  ]= -(d2cdx2[0,:]-cold[0,:]*cold[1,:])
 
    sma[1,:,1]=  ones((1,NJ))
    smd[1,:,0]=  cold[1,:]
    smd[1,:,1]=  cold[0,:]
    smg[1,:  ]= -(d2cdx2[1,:]+cold[0,:]*cold[1,:])
 
    #__________________________________________________  Boundary-Condition 1
    #   sme = smd[:,0,:]
    #   smp = smb[:,0,:]
    #   smf = smg[:,0]
 
    sme = zeros([N,N])
    smp = zeros([N,N])
    smf = zeros([N,1])
 
    sme[0,0] = 1.0
    smf[0  ] = 1.0-cold[0,0]
 
    sme[1,1] = 1.0
    smf[1  ] = 1.0-cold[1,0]
 
    # Insert (sme smp smf) into (smb smd smg)
    smb[:,0,:] = smp[:,:]
    smd[:,0,:] = sme[:,:]
    smg[:,0] = transpose(smf)
 
 
    #_______________________________________________  Boundary-Condition 2
    #   sme = smd[:,NJ-1,:]
    #   smp = smb[:,NJ-1,:]
    #   smf = smg[:,NJ-1]
 
    sme = zeros([N,N])
    smp = zeros([N,N])
    smf = zeros([N,1])
 
    sme[0,0] = 1.0
    smf[0  ] = -cold[0,NJ-1]
 
 
    smp[1,1] = 1.0
    smf[1  ] = -dcdx[1,NJ-1]
 
    # Insert (sme smp smf) into (smb smd smg)
    smb[:,NJ-1,:] = smp[:,:]
    smd[:,NJ-1,:] = sme[:,:]
    smg[:,NJ-1] = transpose(smf)
 
    return ( sma, smb, smd, smg )

def samplefd(DoPlot=True):
    # def samplefd(doplot):
    # ---------------------------------
    #   solves two coupled equations: 
    #    d2cdx2 - c T = 0
    #    d2Tdx2 + c T = 0
    #    with bcs:  at x=0 c=1, T = 1
    #               at x=1 c=0, dTdx=0
    # ________________________________________


    N=2
    NJ=100
    itmax=100
    tol=1e-6
    xmax=1.0

    h=xmax/float(NJ-1)

    cold=zeros((N,NJ))
    delc=cold

    # initial_guess
    cold[0,:]=ones((1,NJ))
    cold[1,:]=zeros((1,NJ))

    cold= bound_val(cold, 
                    deqfun=fillmat,
                    order=2,
                    h=h,
                    tol=tol,
                    itmax=itmax) 

    # print_results;

    x=linspace(0.0,1.0,NJ)
    c=cold[0,:]
    T=cold[1,:]

    if DoPlot:
        import matplotlib.pyplot as mpl
        mpl.plot(x,c,'-o' )
        mpl.plot(x,T,'-s' )
        mpl.legend(['c','T'],loc='center right')
        mpl.xlabel('x / L');
        
    return x, c, T

if __name__ == "__main__":
    samplefd(DoPlot=True)