import numpy as np
import pandas as pd
from pressure import pressure
swi =  0.2 # saturation of water inital 
phi =  0.2 # porocity
tou =   3.7
t =  100000 # time 
dt =  25 # time step 
dx=  0.01 # node
k=  0.94*10**(-12) # permeability 
pini = (1.4*(10**(-7))) # pressure initial
qw= 9*10**(-6) # injection rate
Area= 0.0625 # area of res.
sor =  0.2 # residual oil saturation 
swir =  0.2 # irresidual water saturation
n = 10 # number of nodes
prA = (pini)*np.ones(n,dtype=object) # pressure 
srA = (0.21)*np.ones(n,dtype=object) # saturation 
srAt= (0.21)*np.ones(n,dtype=object) # ittrative step saturation

error=np.zeros(n)
for i in range(0,0.2,0.0137): # time step 
    srAt = srA  # time step saturtion 
    prAt = prA  # 
    pcod= (k*(srAt**(-1/tou))) # old cap pressure 
    while(True):
        srAit=srA  # itrative step saturation 
        pcnew=(k*(srAit**((-1)/tou))) # capillary pressure
        prA,lamw = pressure(srA,prAt,pcnew,pcod,n)  # capillary pressure old time step 
        print("_____________________________________")
        print(prA)
        print("_____________________________________")
        a = prA
        # boundry conditions , first node ,last node , remaining  
        for i in range(n) :
            if i == 0:
                srA[i]=1-sor
            elif i == n-1:
                srA[i]=(srAt[i]-(((lamw[i]+lamw[i-1])*dt/4*phi*(dx**2))*(a[i]+prAt[i]-a[i-1]-prAt[i-1])))
            else:
                srA[i]=((srAt[i]-((((lamw[i]+lamw[i-1])*dt/4*phi*(dx**2))*(a[i]+prAt[i]))-(((lamw[i]+lamw[i+1])*dt/4*phi*(dx**2))*(a[i]+prAt[i]))+(((lamw[i]+lamw[i-1])*dt/4*phi*(dx**2))*(a[i-1]+prAt[i-1]))+(((lamw[i]+lamw[i+1])*dt/4*phi*(dx**2))*(a[i+1]+prAt[i+1])))))
        # to calculate error
        for k in range(n):
            error[k]=srA[k]-srAit[k]
        # to verify 
        if max(error)<(10**(-5)):
            print(srA)
            break  
