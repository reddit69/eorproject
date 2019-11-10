import numpy as np 
pce=3000
muw=0.001
muo=0.000392
qw=9*10**(-6)
Area=0.0625
comp=7*(10**(-8))
con = 3.1731*10**(-7) # change
swi = 0.2
phi = 0.2
tou = 3.7
t = 100.0
dt = 0.025
dx= 0.01
k= 0.94*10**(-12)
pini = (1.4*(10**(-7)))
qw=9*10**(-6)
sor = 0.2
swir = 0.2
n =10
def pressure(srA,prAt,pcnew,pcod,n):
    srAstr = np.ones(n,dtype=object) # saturation 
    coffpnt = np.ones(n,dtype=object) # coff of point 
    kro=np.ones(n,dtype=object) # relative perm. oil 
    krw=np.ones(n,dtype=object) # relative perm. water 
    lamw=np.ones(n,dtype=object) # mobility of water 
    lamo=np.ones(n,dtype=object) # mobility of oil 
    lamt=np.ones(n,dtype=object)  # mobility total
    for i in range(n): # brooks and corey model half implicit , half explicit 
        srAstr[i] = (srA[i]-swir)/(1-sor-swir)
        krw[i]=(srAstr[i]**((2+(3*tou))/(tou)))
        kro[i]=(((1-srAstr[i])**2)*((1-((srAstr[i])**(2+tou)/tou))))
        lamw[i]=((k*krw[i])/muw)
        lamo[i]=((k*kro[i])/muo)
        lamt[i]=lamo[i]+lamw[i] 
        if i == 0 :
            coffpnt[i]=  ((phi*comp*dx/dt)+((lamt[i]+lamt[i+1])/4*(dx)))
        elif i == n-1 :
            coffpnt[i]= ((phi*comp*dx/dt)+((lamt[i]+lamt[i-1])/4*(dx)))
        else :
            coffpnt[i]= ((phi*comp*dx/dt)+(lamt[i]/2*(dx)+((lamt[i-1]+lamt[i+1])/4*(dx))))

    A = np.zeros((n,n),dtype=float) # coff matrix of left side unkonws
    for i in range(n): # calculation of coff  matrix with boundry condition 
        A[i][i] = coffpnt[i]
        if i == 0 :
            A[i][i+1]= -((lamt[i]+lamt[i+1])/(4*(dx)))
        elif i == n-1 :
            A[i][i-1]= -((lamt[i]+lamt[i-1])/(4*(dx)))
        else :
            A[i][i+1]= -((lamt[i]+lamt[i+1])/(4*(dx)))
            A[i][i-1]= -((lamt[i]+lamt[i-1])/(4*(dx)))
    # we are calc left side m
    coffB = np.ones(n) # coff of right side known matrix  // previous time step 
    for i in range(n):
        if i==0 :
            coffB[i] =  (prAt[i]*((phi*comp*dx/dt)-((lamt[i]+lamt[i+1])/4*(dx))))+((lamt[i]+lamt[i+1])*prAt[i+1]/4*(dx)) +(((lamo[i]+lamo[i+1])/4*(dx))*(pcnew[i+1]-pcnew[i]+pcod[i+1]-pcod[i]))-con
        elif i == n-1 :
            coffB[i] = (prAt[i]*((phi*comp*dx/dt)-((lamt[i]+lamt[i-1])/4*(dx)))+((lamt[i]+lamt[i-1])*prAt[i-1]/4*(dx)) -(((lamo[i]+lamo[i-1])/4*(dx))*(pcnew[i]-pcnew[i-1]+pcod[i]-pcod[i-1])))
        else :
            coffB[i] = (((prAt[i]*((phi*comp*dx/dt)-((lamt[i]+lamt[i-1])/4*(dx))-(lamt[i]+lamt[i+1])/4*(dx)))+((lamt[i]+lamt[i-1])*prAt[i-1]/4*(dx))+((lamt[i]+lamt[i+1])*prAt[i+1]/4*(dx)) +(((lamo[i]+lamo[i-1])/4*(dx))-((lamo[i]+lamo[i+1])/4*(dx)))*(pcnew[i]+pcod[i]))+((pcnew[i+1]+pcod[i+1])*((lamo[i]+lamo[i+1])/4*(dx)))+((pcnew[i-1]+pcod[i-1])*((lamo[i]+lamo[i-1])/4*(dx))))

    a = np.linalg.inv(A) # 
    x = np.matmul(a,coffB)
    return x.transpose(), lamw 
