import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd

#Hestro nov. 2020 (adapted from Dullemond uni. Heidelberg)
#Modification may 2021 by PAJUELO Pierre

#------------------------------------------------------
#Définitions des fonctions
#PFD
def f(y,t,mu):
    x = y[0:3].copy()
    v = y[3:6].copy()
    r = np.sqrt(x[0]**2+x[1]**2+x[2]**2)
    dxdt = v
    dvdt = (-mu*x)/(r**3)
    dy = np.hstack((dxdt,dvdt))
    return dy
#Norme d'un vecteur
def norm(x):
    a,b,c=x
    z=np.sqrt(a**2+b**2+c**2)
    return z

#Valeur moyenne d'une liste
def moy(A):
    if len(A)==0:
        return 0
    x=A[0]
    for i in range(1,len(A),1):
        x=x+A[i]
    return x/len(A)
#-----------------------------------------------
# Central force
mstar = 2*10**30
mu = 6.67*10**(-11)*mstar                  # Gravitational parameter
au = 149597870.7*10**3 #meters

# Time values for orbit propagation
jour = 86400                # Time Value
to = 0
#t0 = 4*365*jour
#t0 = (2.66/2)*365*jour                      # Starting time
a = ((0.89+2.94)/2)*au
t0 = np.sqrt((a**3*4*np.pi**2)/mu)/2 #3eme loi de Kepler pour obtenir la demi période
ti = 20*jour               # Impulsion frequency
tendi = t0+ 1000*ti              # End of Impulsion
na = int((tendi-t0)/ti)          # Echantillonage Impulsion
nb = int(t0/ti)
tend= tendi + 10*365*jour            # End time
nt = 100                   # Echantillonage Integration
t = np.linspace(t0,ti,nt)   # Vector list for values of time
# Example : tfinal=np.linspace(t0,tend,nt) 
#-----------------------------------------------
#Sommaire
protocole1='Désactivé'
etude1='Non'
etude2='Non'
OrbA='Non'
OrbB='Non'
OrbC='Non'
protocole2='Désactivé'
protocole3='Désactivé'
etude3='Non'
etude4='Non'
OrbD='Non'
OrbE='Non'
OrbF='Non'
protocole4='Désactivé'
protocole5='Désactivé'
base='Non'
etude5='Non'
etude6='Non'
etude7='Non'
OrbG='Non'
OrbH='Non'
OrbI='Non'
protocole6='Désactivé'
protocole7='Désactivé'
protocole8='Désactivé'
asteroidonly='Non'
solarsystem='Non'
protocoleTerre='Désactivé'
plot='Non'
#-----------------------------------------------
if protocole1=='Activé':
    # Protocole 1 : Impulsion discrétisée
    
    #Initial conditions at time t0
    r0 = 2.94*au
    #v0 = 0.5*29.5*10**3 #m.s-1
    e = (2.94-0.89)/(2.94+0.89)
    a = ((0.89+2.94)/2)*au
    v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
    Tasteroide = np.sqrt((a**3*4*np.pi**2)/mu)
    #Example : np.sqrt(mu/r0) or 29.5*10**3
    #ODG : km.s-1=v0 to mm.s-1=deltav and tend=qqs dizaines années    
    x0 = np.array([r0,0.,0.],dtype='float64') # initial condition position
    v0 = np.array([0.,v0,0.])                 # initial condition velocity
    y0 = np.hstack((x0,v0))                   # initial condition state-vector
    
    #Impulsion parameters 
    #deltav = 5*10**(0)                        # Value of Impulsion 
    #u = np.array([1.,1.,0.],dtype='float64') # Impulsion vector
    
    #Relation de récurrence de la vitesse : vi+1 = vi + deltav(cte)*vi(vecteur unitaire)
    
    #Définition des constantes de temps
    to=0
    t0=Tasteroide
    print(t0/jour)
    nb=200
    impulsion=30*jour
    tiimpulsion=t0-impulsion/2
    tfinal1=np.linspace(to,tiimpulsion,nb)       # Timeline impulsion
    #ti=0.3*jour
    TFINAL1=tiimpulsion-to
    pas=TFINAL1/nb
    
    tendi=t0+impulsion/2
    TFINAL2=tendi-(tiimpulsion+pas)
    na=int(TFINAL2/pas)
    tfinal2=np.linspace(tiimpulsion+pas,tendi,na)
    
    tend=tendi+50*t0
    TFINAL3=tend-(tendi+pas)
    nc=int(TFINAL3/pas)
    tfinal3=np.linspace(tendi+pas,tend,nc)
    tfinal4=np.concatenate((tfinal1,tfinal2,tfinal3))
    n = len(tfinal4)
    nt=200
    j=0
    
    Etude=tendi+4*t0
    
    Rf=[[],[],[],[]]
    GOAL=2
    DeltaV=np.linspace(0,100,GOAL)
    Compteur=0
    Compteurgoal=2
    Compteuri=1
    Compteurgoal1=-1
    Compteuri1=0
    for deltav in DeltaV:
        Rf[0].append(deltav)
        R=[]
        r0 = 2.94*au
    
        e = (2.94-0.89)/(2.94+0.89)
        a = ((0.89+2.94)/2)*au
        v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
        x0 = np.array([r0,0.,0.],dtype='float64') # initial condition position
        v0 = np.array([0.,v0,0.])                 # initial condition velocity
        y0 = np.hstack((x0,v0))                   # initial condition state-vector
    
        for i in tfinal4:
            if i==to: #échelle de temps de l'impulsion
                
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                v0 = sol[:,3:6].T
                if OrbA=='Oui':
                    if Compteuri<=Compteur and Compteur<=Compteurgoal:
                        orb1=plt.plot(x0[0,:],x0[1,:],'g',label='Orbite avant impulsion discrétisée')
                        j+=1
                        #print(1)
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                y0 = np.hstack((x1,v1))
                
            if i<tiimpulsion:
                
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                
                v0 = sol[:,3:6].T
                if OrbA=='Oui':
                    if Compteuri<=Compteur and Compteur<=Compteurgoal:
                        orb1=plt.plot(x0[0,:],x0[1,:],'g')
                        j+=1
                        #print(2)
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                y0 = np.hstack((x1,v1))
                
            if i==tiimpulsion: #échelle de temps de l'impulsion
                
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                v0 = sol[:,3:6].T
                if OrbB=='Oui':
                    if Compteuri<=Compteur and Compteur<=Compteurgoal:
                        orb1=plt.plot(x0[0,:],x0[1,:],'r',label='Impulsion discrétisée')
                        j+=1
                        #print(3)
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                r=norm(x1)
                
                
                a,b,c=x1/r
                utheta=np.array([-b,a,0.])
                v=norm(v1)
                u=np.vdot(v1,utheta)*utheta/v
                #u=v1/v
                #u=x1/r
                #print(u)
                #print(x1)
                #print(r)
                v1 = np.add(v1,deltav*u)                             # IMPULSION
                y0 = np.hstack((x1,v1))
                
            if i>tiimpulsion and i<=tendi: #échelle de temps de l'impulsion
                
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                v0 = sol[:,3:6].T
                if OrbB=='Oui':
                    if Compteuri<=Compteur and Compteur<=Compteurgoal:
                        orb1=plt.plot(x0[0,:],x0[1,:],'r')
                        j+=1
                        #print(4)
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                r=norm(x1)
                a,b,c=x1/r
                utheta=np.array([-b,a,0.])
                v=norm(v1)
                u=np.vdot(v1,utheta)*utheta/v
                #print(utest)
                #print(deltav)
                #u=x1/r
                v1 = np.add(v1,deltav*u)                             # IMPULSION
                y0 = np.hstack((x1,v1))
                
            if i>tendi:
                
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                v0 = sol[:,3:6].T
                if OrbC=='Oui':
                    if Compteuri<=Compteur and Compteur<=Compteurgoal:
                        orb1=plt.plot(x0[0,:],x0[1,:],'b')
                        j+=1
                        #print(5)
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                y0 = np.hstack((x1,v1)) 
                r=norm(x1)
                if r not in R:
                    if i>Etude:
                        R.append(r)
                
            if  i==tfinal4[n-1]:
                
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                v0 = sol[:,3:6].T
                if OrbC=='Oui':
                    if Compteuri<=Compteur and Compteur<=Compteurgoal:
                        orb1=plt.plot(x0[0,:],x0[1,:],'b',label='Après impulsion discrétisée')
                        j+=1
                        #print(6)
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                y0 = np.hstack((x1,v1)) 
                r=norm(x1)
                if r not in R:
                    if i>Etude:
                        R.append(r)
        if Compteuri1<=Compteur and Compteur<=Compteurgoal1:
            Rlen=len(R)
            Temps=[i for i in range(Rlen)]
            plt.plot(Temps,R)
            
        Compteur+=1
        if (Compteur/GOAL*100)%10==0:
            print('Etat :',Compteur/GOAL*100)
        Rmin=min(R)
        Rmax=max(R)
        Rf[1].append(Rmin)
        Rf[2].append(Rmax)
        a=(Rmax+Rmin)/2
        Rf[3].append(np.sqrt((a**3*4*np.pi**2)/mu))
    #plt.show()
    if etude1=='Oui':    
        plt.plot(Rf[0],Rf[1])
        U = {'Delta v': Rf[0],'Raphelie': Rf[1]}
        df = pd.DataFrame(U, columns = ['Delta v','Raphelie'])
        df.to_excel (r'C:\Users\pierr\Downloads\Regression.xlsx', index = False,columns=['Delta v', 'Raphelie'])
        plt.xlabel("Deltav")
    
        plt.ylabel("Rpérihélie")
        plt.show()
    if etude2=='Oui':
        plt.plot(Rf[0],Rf[3])
        plt.xlabel("Deltav")
        plt.ylabel("Période")
        U = {'Delta v': Rf[0],'Période': Rf[3]}
        df = pd.DataFrame(U, columns = ['Delta v','Période'])
        df.to_excel (r'C:\Users\pierr\Downloads\Regression.xlsx', index = False,columns=['Delta v', 'Période'])
        plt.show()
        
    
    

#-----------------------------------------------
if protocole2=='Activé':
    # Protocole 2 : Orbite de base 
    
    #Conditions initiales
    r0 = 2.94*au
    #v0 = 1.3*np.sqrt(mu/r0)
    a = ((0.89+2.94)/2)*au
    e = (2.94-0.89)/(2.94+0.89)
    v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
    x0 = np.array([r0,0.,0.],dtype='float64')
    v0 = np.array([0.,v0,0.],dtype='float64')
    y0 = np.hstack((x0,v0))
    #Orbite de base
    t0=0
    tend=3*365*jour
    nt=500000
    tfinal=np.linspace(t0,tend,nt)
    sol = odeint(f,y0,tfinal,args=(mu,))
    x0 = sol[:,0:3].T
    v0 = sol[:,3:6].T
    orb2=plt.plot(x0[0,:],x0[1,:],'m',label='Orbite de départ')

#-----------------------------------------------

if protocole3=='Activé':

    # Protocole 3 : Impulsion singulière

    
    #Conditions initiales
    deltav = 1*10**(3)
    r0 = 2.94*au
    a = ((0.89+2.94)/2)*au
    e = (2.94-0.89)/(2.94+0.89)
    v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
    x0 = np.array([r0,0.,0.],dtype='float64')
    v0 = np.array([0.,v0,0.],dtype='float64')
    y0 = np.hstack((x0,v0))
    Tasteroide = np.sqrt((a**3*4*np.pi**2)/mu)
    #Structuration temporelle
    to=0
    t0=Tasteroide
    nb=1500
    

    tfinal1=np.linspace(to,t0,nb)       # Timeline impulsion
    TFINAL1=t0-to
    pas=TFINAL1/nb
    
    tend=t0+5*t0
    Etude=3*t0
    TFINAL2=tend-(t0+pas)
    na=int(TFINAL2/pas)
    
    tfinal2=np.linspace(t0+pas,tend,na)
    tfinal4=np.concatenate((tfinal1,tfinal2))
    n = len(tfinal4)
    
    
    #Constantes
    nt=200
    GOAL=100
    DELTAV=np.linspace(0,10**5,GOAL)
    Rf=[[],[],[],[]]
    Compteur=0
    #Impulsion distincte (pendant une durée ti)
    for deltav in DELTAV:
        R=[]
        Rf[0].append(deltav)
        for i in tfinal4:
            if i==to: #échelle de temps de l'impulsion
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                v0 = sol[:,3:6].T
                if OrbD=='Oui':
                    orb1=plt.plot(x0[0,:],x0[1,:],'g',label='Orbite avant impulsion distincte')
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
        
                y0 = np.hstack((x1,v1))
                
            if i<t0:
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                
                v0 = sol[:,3:6].T
                if OrbD=='Oui':
                    orb1=plt.plot(x0[0,:],x0[1,:],'g')
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                y0 = np.hstack((x1,v1))
                   
            if i==t0: #échelle de temps de l'impulsion
                t = np.linspace(t0,t0+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                v0 = sol[:,3:6].T
                if OrbE=='Oui':
                    orb3=plt.plot(x0[0,:],x0[1,:],'y',label='Impulsion distincte')
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                
                r=norm(x1)
                a,b,c=x1/r
                utheta=np.array([-b,a,0.])
                v=norm(v1)
                u=np.vdot(v1,utheta)*utheta/v
                v1 = np.add(v1,deltav*u)
                y0 = np.hstack((x1,v1))
            if i>t0 and i!=tfinal4[n-1]:
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                v0 = sol[:,3:6].T
                if OrbF=='Oui':
                    orb3=plt.plot(x0[0,:],x0[1,:],'b')
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                y0 = np.hstack((x1,v1))
                r=norm(x1)
                if r not in R:
                        if i>Etude:
                            R.append(r)
            if i==tfinal4[n-1]:
                t = np.linspace(i,i+pas,nt)
                sol = odeint(f,y0,t,args=(mu,))
                x0 = sol[:,0:3].T
                v0 = sol[:,3:6].T
                if OrbF=='Oui':
                    orb3=plt.plot(x0[0,:],x0[1,:],'b',label='Courbe après impulsion distincte')
                x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                y0 = np.hstack((x1,v1))
        Rmin=min(R)
        Rf[1].append(Rmin)
        Compteur+=1
        if (Compteur/GOAL*100)%10==0:
            print('Etat :',Compteur/GOAL*100)
        
        Rmax=max(R)
        
        Rf[2].append(Rmax)
        a=(Rmax+Rmin)/2
        Rf[3].append(np.sqrt((a**3*4*np.pi**2)/mu))
    if etude3=='Oui':    
        plt.plot(Rf[0],Rf[1])
        U = {'Delta v': Rf[0],'Raphelie': Rf[1]}
        df = pd.DataFrame(U, columns = ['Delta v','Raphelie'])
        df.to_excel (r'C:\Users\pierr\Downloads\Regression.xlsx', index = False,columns=['Delta v', 'Raphelie'])
        plt.xlabel("Deltav")
        plt.ylabel("Rpérihélie")
        plt.show()
    if etude4=='Oui':
        plt.plot(Rf[0],Rf[3])
        plt.xlabel("Deltav")
        plt.ylabel("Période")
        U = {'Delta v': Rf[0],'Période': Rf[3]}
        df = pd.DataFrame(U, columns = ['Delta v','Période'])
        df.to_excel (r'C:\Users\pierr\Downloads\Regression.xlsx', index = False,columns=['Delta v', 'Période'])
        plt.show()
#-----------------------------------------------


if protocole4=='Activé':
    #Protocole 4 : Calcul des points de collision
    np.seterr(divide='ignore', invalid='ignore')
    #Conditions initiales
    r1=0.89
    r2=2.94
    #Astéroïde
    r0 = max(r1,r2)*au
    a = ((r1+r2)/2)*au
    e = (abs(r1-r2))/(r1+r2)
    v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
    
    
    #Initialisation
    Tasteroide = np.sqrt((a**3*4*np.pi**2)/mu)
    x0 = np.array([r0,0.,0.],dtype='float64')
    v0 = np.array([0.,v0,0.],dtype='float64')
    y0 = np.hstack((x0,v0))
    
    #Terre
    r = 152*10**9
    v1 = 29.78*10**3
    Tterre = np.sqrt((r**3*4*np.pi**2)/mu)

    #Initialisation
    x2 = np.array([r,0.,0.],dtype='float64')
    v2 = np.array([0.,v1,0.],dtype='float64')
    y1 = np.hstack((x2,v2))
    
    #Temps
    to=0
    t0=Tasteroide
    
    nb=200
    impulsion=30*jour
    tiimpulsion=t0-impulsion/2
    tfinal1=np.linspace(to,tiimpulsion,nb)       # Timeline impulsion
    #ti=0.3*jour
    TFINAL1=tiimpulsion-to
    pas=TFINAL1/nb
    
    tendi=t0+impulsion/2
    TFINAL2=tendi-(tiimpulsion+pas)
    na=int(TFINAL2/pas)
    tfinal2=np.linspace(tiimpulsion+pas,tendi,na)
    
    tend=tendi+50*t0
    TFINAL3=tend-(tendi+pas)
    nc=int(TFINAL3/pas)
    tfinal3=np.linspace(tendi+pas,tend,nc)
    tfinal4=np.concatenate((tfinal1,tfinal2,tfinal3))
    n = len(tfinal4)
    
    nt=200
    
    
    TERRE=[[],[],[]]
    ASTEROIDE = [[],[],[]]
    for i in tfinal4:
        if i<=tiimpulsion:
            ti = pas
            t = np.linspace(i,i+ti,nt)
            #Position Terre
            (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
            x1 = sol1[:,0:3].T
            v1 = sol1[:,3:6].T
            plt.plot(x1[0,:],x1[1,:],'b')
            
                
            x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
            v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
            y1 = np.hstack((x2,v2))
            #Position Astéroide
            sol = odeint(f,y0,t,args=(mu,))
            x0 = sol[:,0:3].T
            v0 = sol[:,3:6].T
            orb1=plt.plot(x0[0,:],x0[1,:],'g')
            x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
            v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
            y0 = np.hstack((x1,v1))
            
        if i>tiimpulsion and i<=tendi: #échelle de temps de l'impulsion
            ti = pas
            t = np.linspace(i,i+ti,nt)
            
            
           
           #Position Terre
            (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
            x1 = sol1[:,0:3].T
            v1 = sol1[:,3:6].T
            plt.plot(x1[0,:],x1[1,:],'b')
            
                
            x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
            v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
            y1 = np.hstack((x2,v2))


            sol2 = odeint(f,y0,t,args=(mu,))
            x0 = sol2[:,0:3].T
            v0 = sol2[:,3:6].T
            plt.plot(x0[0,:],x0[1,:],'g')

            x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
            v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
            
            y0 = np.hstack((x1,v1))
        if i>tendi:
            ti = pas
            t = np.linspace(i,i+ti,nt)
            
            
            #Position Terre
            (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
            x1 = sol1[:,0:3].T
            v1 = sol1[:,3:6].T
            plt.plot(x1[0,:],x1[1,:],'b')
            if i<(tendi+Tterre):
                X1 = x1.tolist()
                TERRE[0].extend(X1[0])
                TERRE[1].extend(X1[1])
                TERRE[2].extend(X1[2])
                
            x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
            v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
            y1 = np.hstack((x2,v2))


            sol2 = odeint(f,y0,t,args=(mu,))
            x0 = sol2[:,0:3].T
            v0 = sol2[:,3:6].T
            plt.plot(x0[0,:],x0[1,:],'g')

            if i<(tendi+Tasteroide):
                X0 = x0.tolist()
                ASTEROIDE[0].extend(X0[0])
                ASTEROIDE[1].extend(X0[1])
                ASTEROIDE[2].extend(X0[2])

            x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
            v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
            y0 = np.hstack((x1,v1))
            
    #Eviter les doublons
    TERRE2 = [[],[],[]] 
    for i in TERRE[0] : 
        if i not in TERRE2[0]: 
            TERRE2[0].append(i) 
    for i in TERRE[1] : 
        if i not in TERRE2[1]: 
            TERRE2[1].append(i) 
    for i in TERRE[2] : 
        if i not in TERRE2[2]: 
            TERRE2[2].append(i) 
    ASTEROIDE2 = [[],[],[]] 
    for i in ASTEROIDE[0] : 
        if i not in ASTEROIDE2[0]: 
            ASTEROIDE2[0].append(i) 
    for i in ASTEROIDE[1] : 
        if i not in ASTEROIDE2[1]: 
            ASTEROIDE2[1].append(i) 
    for i in ASTEROIDE[2] : 
        if i not in ASTEROIDE2[2]: 
            ASTEROIDE2[2].append(i)
            

    print('Nombre de points :', 'Terre :',len(TERRE2[0]),'Astéroide :',len(ASTEROIDE2[0]))
    for x in TERRE2[0]:
        for y in ASTEROIDE2[0]:
            if abs(x-y)<10**(7):
                #print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                p=TERRE2[1][TERRE2[0].index(x)]
                c=ASTEROIDE2[1][ASTEROIDE2[0].index(y)]
                #print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                
                if abs(p-c)<10**(7):
                    
                    print('CRASH !!')
                    print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                    print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                    
                    plt.plot(x,p,'o')
                    plt.plot(y,c,'o',label='Collision')
#-----------------------------------------------
import time


if protocole5=='Activé':
    if etude6=='Oui':
        #Etude 6 : Etude d'optimisation temporelle
        
        discretisationf=np.linspace(10,100,20)
        ETUDE6=[[],[]]
        for discretisation in discretisationf:
            ETUDE6[0].append(int(discretisation))
            print('Statut :',int(discretisation),'sur',100)
            start_time=time.time()
            #Protocole 5 : Modification de la trajectoire de l'astéroïde
            np.seterr(divide='ignore', invalid='ignore')
            #Conditions initiales
            
            r1=0.89
            r2=2.94
            #Astéroïde
            
            r0 = max(r1,r2)*au
            a = ((r1+r2)/2)*au
            e = (abs(r1-r2))/(r1+r2)
            v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
            #print(v0)
            
            #Initialisation
            Tasteroide = np.sqrt((a**3*4*np.pi**2)/mu)
            x0 = np.array([r0,0.,0.],dtype='float64')
            v0 = np.array([0.,v0,0.],dtype='float64')
            y0 = np.hstack((x0,v0))
            
            #Terre
            r = 152*10**9
            v1 = 29.78*10**3
            Tterre = np.sqrt((r**3*4*np.pi**2)/mu)
        
            #Initialisation
            x2 = np.array([r,0.,0.],dtype='float64')
            v2 = np.array([0.,v1,0.],dtype='float64')
            y1 = np.hstack((x2,v2))
        
            #Echelle de Période
            #PERIODE=np.linspace(1,periodemax,nbp)
            PARAMETRE=[[],[]]
            PARAMETRE2=[[],[]]
            #for periode in PERIODE:
            #print('Periode actuel :',periode)
            #Temps
            to=0
            t0=Tasteroide
            
            nb=int(discretisation)
            impulsion=30*jour
            tiimpulsion=t0-impulsion/2
            tfinal1=np.linspace(to,tiimpulsion,nb)       # Timeline impulsion
            
            TFINAL1=tiimpulsion-to
            pas=TFINAL1/nb
            
            tendi=t0+impulsion/2
            TFINAL2=tendi-(tiimpulsion+pas)
            na=int(TFINAL2/pas)
            tfinal2=np.linspace(tiimpulsion+pas,tendi,na)
            
            tend=tendi+2*t0
            TFINAL3=tend-(tendi+pas)
            nc=int(TFINAL3/pas)
            tfinal3=np.linspace(tendi+pas,tend,nc)
            tfinal4=np.concatenate((tfinal1,tfinal2,tfinal3))
            n = len(tfinal4)
            #Données du protocole:
            nt=200
            critere=8
            critereetude=9
            
            deltavmax=1000
            precision=2
            vfinal = np.linspace(0.,deltavmax,precision)
            
            
            U=[[],[]]
            P=[[0 for i in range(precision)] for j in range(5)]
            P2=[[0 for i in range(precision)] for j in range(5)]
            
            indice=0
            for deltav in vfinal:
                TERRE=[[],[],[]]
                ASTEROIDE = [[],[],[]]
                ASTEROIDE3 = [[],[],[]]
                A=[[],[]]
                A2=[[],[]]
                for i in tfinal4:
                    if i<tiimpulsion:
                        ti = pas
                        t = np.linspace(i,i+ti,nt)
                        
                        #Position Terre
                        (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
                        x1 = sol1[:,0:3].T
                        v1 = sol1[:,3:6].T
                        if OrbG=='Oui':
                            plt.plot(x1[0,:],x1[1,:],'b')
                        
                        x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
                        v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
                        y1 = np.hstack((x2,v2))
                        
                        #Position Astéroide
                        sol = odeint(f,y0,t,args=(mu,))
                        x0 = sol[:,0:3].T
                        v0 = sol[:,3:6].T
                        if OrbG=='Oui':
                            orb1=plt.plot(x0[0,:],x0[1,:],'g')
                        x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                        v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                        y0 = np.hstack((x1,v1))
                
                        
                    if i>=tiimpulsion and i<=tendi: #échelle de temps de l'impulsion
                        ti = pas
                        t = np.linspace(i,i+ti,nt)
                        
                        #Position Terre
                        (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
                        x1 = sol1[:,0:3].T
                        v1 = sol1[:,3:6].T
                        if OrbH=='Oui':
                            plt.plot(x1[0,:],x1[1,:],'b')
                        
                        x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
                        v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
                        y1 = np.hstack((x2,v2))
                        
                        #Position Astéroide
                        sol = odeint(f,y0,t,args=(mu,))
                        x0 = sol[:,0:3].T
                        v0 = sol[:,3:6].T
                        if OrbH=='Oui':
                            orb1=plt.plot(x0[0,:],x0[1,:],'g')
                        
                        x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                        v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                        #-------------------
                        r=norm(x1)
                        a,b,c=x1/r
                        utheta=np.array([-b,a,0.])
                        '''
                        v=norm(v1)
                        u=np.vdot(v1,utheta)*utheta/v
                        '''
                        #--------------------------
                        v1 = np.add(v1,deltav*utheta)
                        y0 = np.hstack((x1,v1))
                    if i>tendi:
                        ti = pas
                        t = np.linspace(i,i+ti,nt)
                        
                        
                        #Position Terre
                        (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
                        x1 = sol1[:,0:3].T
                        v1 = sol1[:,3:6].T
                        if OrbI=='Oui':
                            plt.plot(x1[0,:],x1[1,:],'b')
                        if i<(tendi+Tterre):
                            X1 = x1.tolist()
                            TERRE[0].extend(X1[0])
                            TERRE[1].extend(X1[1])
                            TERRE[2].extend(X1[2])
                            
                        x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
                        v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
                        y1 = np.hstack((x2,v2))
            
            
                        sol2 = odeint(f,y0,t,args=(mu,))
                        x0 = sol2[:,0:3].T
                        v0 = sol2[:,3:6].T
                        if OrbI=='Oui':
                            plt.plot(x0[0,:],x0[1,:],'g')
            
                        if i<(tendi+Tasteroide):
                            X0 = x0.tolist()
                            ASTEROIDE[0].extend(X0[0])
                            ASTEROIDE[1].extend(X0[1])
                            ASTEROIDE[2].extend(X0[2])
            
                        x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                        v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                        y0 = np.hstack((x1,v1))
        
        
                        if i<(tendi+Tasteroide):
                            X0 = x0.tolist()
                            ASTEROIDE[0].extend(X0[0])
                            ASTEROIDE[1].extend(X0[1])
                            ASTEROIDE[2].extend(X0[2])
                        if etude5=='Oui':
                            periode=5
                            if i>(tendi+(periode-1)*Tasteroide) and i<(tendi+periode*Tasteroide):
                                X0 = x0.tolist()
                                ASTEROIDE3[0].extend(X0[0])
                                ASTEROIDE3[1].extend(X0[1])
                                ASTEROIDE3[2].extend(X0[2])
                            
                        x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                        v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                        y0 = np.hstack((x1,v1))
                #Eviter les doublons
                TERRE2 = [[],[],[]] 
                for i in TERRE[0] : 
                    if i not in TERRE2[0]: 
                        TERRE2[0].append(i) 
                        
                for i in TERRE[1] : 
                    if i not in TERRE2[1]: 
                        TERRE2[1].append(i) 
                for i in TERRE[2] : 
                    if i not in TERRE2[2]: 
                        TERRE2[2].append(i) 
                ASTEROIDE2 = [[],[],[]] 
                for i in ASTEROIDE[0] : 
                    if i not in ASTEROIDE2[0]: 
                        ASTEROIDE2[0].append(i) 
                for i in ASTEROIDE[1] : 
                    if i not in ASTEROIDE2[1]: 
                        ASTEROIDE2[1].append(i) 
                for i in ASTEROIDE[2] : 
                    if i not in ASTEROIDE2[2]: 
                        ASTEROIDE2[2].append(i)
                if etude5=='Oui':
                    ASTEROIDE4 = [[],[],[]] 
                    for i in ASTEROIDE3[0] : 
                        if i not in ASTEROIDE4[0]: 
                            ASTEROIDE4[0].append(i) 
                    for i in ASTEROIDE3[1] : 
                        if i not in ASTEROIDE4[1]: 
                            ASTEROIDE4[1].append(i) 
                    for i in ASTEROIDE3[2] : 
                        if i not in ASTEROIDE4[2]: 
                            ASTEROIDE4[2].append(i)
                parametre=0
                
            
                print('Nombre de points :','Terre:', len(TERRE2[0]),'Astéroide:',len(ASTEROIDE2[0]))
                if etude5=='Oui':
                    print('Nombre de points Etude 5:','Nombre de périodes:',periode,'Terre:',len(TERRE2[0]),'Astéroide:',len(ASTEROIDE4[0]))
                for x in TERRE2[0]:
                    for y in ASTEROIDE2[0]:
                        if abs(x-y)<10**(critere):
                            #print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                            p=TERRE2[1][TERRE2[0].index(x)]
                            c=ASTEROIDE2[1][ASTEROIDE2[0].index(y)]
                            #print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                            
                            if abs(p-c)<10**(critere):
                                '''
                                print('Test actuel :',deltav)
                                print('CRASH !!')
                                print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                                print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                                '''
                                ''' #Principe : Prendre la première valeur
                                if x<0 and P[2][indice]==0:
                                    P[0][indice]=x
                                    P[4][indice]=deltav
                                if p>0 and P[3][indice]==0:
                                    P[1][indice]=p
                                '''
                                if x<0:
                                    A[0].append(x)
                                    
                                if p>0:
                                    A[1].append(p)
                                #plt.plot(x,p,'o')
                                #plt.plot(y,c,'o')
                                parametre+=1
                if etude5=='Oui':                
                    parametre2=0
                    for x in TERRE2[0]:
                        for y in ASTEROIDE4[0]:
                            if abs(x-y)<10**(critereetude):
                                #print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                                p=TERRE2[1][TERRE2[0].index(x)]
                                c=ASTEROIDE4[1][ASTEROIDE4[0].index(y)]
                                #print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                                
                                if abs(p-c)<10**(critereetude):
                                    
                                    if x<0:
                                        A2[0].append(x)
                                        
                                    if p>0:
                                        A2[1].append(p)
                                    #plt.plot(x,p,'o')
                                    #plt.plot(y,c,'o')
                                    parametre2+=1
                P[0][indice]=moy(A[0])
                P[1][indice]=moy(A[1])
                P[4][indice]=deltav
                
                if etude5=='Oui':
                    P2[0][indice]=moy(A2[0])
                    P2[1][indice]=moy(A2[1])
                    P2[4][indice]=deltav
                indice+=1
                print('Etat actuel:',indice,'Goal',precision)
                if parametre==0 and deltav!=0:
                    print('Test actuel :',deltav,'Pas de crash')
                    PARAMETRE[0].append(parametre)
                    PARAMETRE[1].append(deltav)
                    
                if etude5=='Oui':
                    if parametre2==0 and deltav!=0:
                        print('Test actuel après plusieurs périodes:',deltav,'Pas de crash')
                        PARAMETRE2[0].append(parametre2)
                        PARAMETRE2[1].append(deltav)
                U[0].append(parametre)
                U[1].append(deltav)
                
         
            #plt.plot(TERRE2[0],TERRE2[1],label='Orbite de la Terre')
            #plt.plot(P[0],P[1],label='Déplacement du crash 1ere periode')
            if etude5=='Oui':
                plt.plot(P2[0],P2[1],label='Déplacement du crash derniere periode')
            
            print('1ère période',PARAMETRE[0][0],PARAMETRE[1][0])
            if etude5=='Oui':
                print('Dernière période',PARAMETRE2[0][0],PARAMETRE2[1][0])
            ETUDE6[1].append(time.time()-start_time)
        plt.plot(ETUDE6[0],ETUDE6[1])
        U = {'Discretisation': ETUDE6[0],'Temps': ETUDE6[1]}
        df = pd.DataFrame(U, columns = ['Discretisation','Temps'])
        df.to_excel(r'C:\Users\pierr\Downloads\Regression.xlsx', index = False,columns=['Discretisation', 'Temps'])
    
    if base=='Oui':
        #Protocole 5 : Modification de la trajectoire de l'astéroïde
        np.seterr(divide='ignore', invalid='ignore')
        #Conditions initiales
        
        r1=0.89
        r2=2.94
        #Astéroïde
        
        r0 = max(r1,r2)*au
        a = ((r1+r2)/2)*au
        e = (abs(r1-r2))/(r1+r2)
        v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
        #print(v0)
        
        #Initialisation
        Tasteroide = np.sqrt((a**3*4*np.pi**2)/mu)
        x0 = np.array([r0,0.,0.],dtype='float64')
        v0 = np.array([0.,v0,0.],dtype='float64')
        y0 = np.hstack((x0,v0))
        
        #Terre
        r = 152*10**9
        v1 = 29.78*10**3
        Tterre = np.sqrt((r**3*4*np.pi**2)/mu)
        
        #Initialisation
        x2 = np.array([r,0.,0.],dtype='float64')
        v2 = np.array([0.,v1,0.],dtype='float64')
        y1 = np.hstack((x2,v2))
        
        #Echelle de Période
        #PERIODE=np.linspace(1,periodemax,nbp)
        PARAMETRE=[[],[]]
        PARAMETRE2=[[],[]]
        #for periode in PERIODE:
        #print('Periode actuel :',periode)
        #Temps
        to=0
        t0=Tasteroide
        
        nb = 1000
        impulsion=30*jour
        tiimpulsion=t0-impulsion/2
        tfinal1=np.linspace(to,tiimpulsion,nb)       # Timeline impulsion
        
        TFINAL1=tiimpulsion-to
        pas=TFINAL1/nb
        
        tendi=t0+impulsion/2
        TFINAL2=tendi-(tiimpulsion+pas)
        '''
        na=int(TFINAL2/pas)
        print('na :', na)
        '''
        na=int(TFINAL2/(0.1*jour))
        print("na:",na)
        tfinal2=np.linspace(tiimpulsion+pas,tendi,na)
        
        tend=tendi+2*t0
        TFINAL3=tend-(tendi+pas)
        '''
        nc=int(TFINAL3/pas)
        print('nc :', nc)
        '''
        nc = 40000
        tfinal3=np.linspace(tendi+pas,tend,nc)
        tfinal4=np.concatenate((tfinal1,tfinal2,tfinal3))
        n = len(tfinal4)
        #Données du protocole:
        nt=2
        critere=8
        critereetude=9
        
        deltavmax=2.3
        precision=5
        vfinal = np.linspace(2,deltavmax,precision)
        
        
        U=[[],[]]
        P=[[0 for i in range(precision)] for j in range(5)]
        P2=[[0 for i in range(precision)] for j in range(5)]
        
        indice=0
        for deltav in vfinal:
            TERRE=[[],[],[]]
            ASTEROIDE = [[],[],[]]
            ASTEROIDE3 = [[],[],[]]
            A=[[],[]]
            A2=[[],[]]
            for i in tfinal4:
                if i<tiimpulsion:
                    ti = TFINAL1/nb
                    t = np.linspace(i,i+ti,nt)
                    
                    #Position Terre
                    (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
                    x1 = sol1[:,0:3].T
                    v1 = sol1[:,3:6].T
                    if OrbG=='Oui':
                        plt.plot(x1[0,:],x1[1,:],'b')
                    
                    x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
                    v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
                    y1 = np.hstack((x2,v2))
                    
                    #Position Astéroide
                    sol = odeint(f,y0,t,args=(mu,))
                    x0 = sol[:,0:3].T
                    v0 = sol[:,3:6].T
                    if OrbG=='Oui':
                        orb1=plt.plot(x0[0,:],x0[1,:],'g')
                    x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                    v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                    y0 = np.hstack((x1,v1))
            
                    
                if i>=tiimpulsion and i<=tendi: #échelle de temps de l'impulsion
                    ti = TFINAL2/na
                    t = np.linspace(i,i+ti,nt)
                    
                    #Position Terre
                    (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
                    x1 = sol1[:,0:3].T
                    v1 = sol1[:,3:6].T
                    if OrbH=='Oui':
                        plt.plot(x1[0,:],x1[1,:],'b')
                    
                    x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
                    v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
                    y1 = np.hstack((x2,v2))
                    
                    #Position Astéroide
                    sol = odeint(f,y0,t,args=(mu,))
                    x0 = sol[:,0:3].T
                    v0 = sol[:,3:6].T
                    if OrbH=='Oui':
                        orb1=plt.plot(x0[0,:],x0[1,:],'g')
                    
                    x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                    v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                    #-------------------
                    r=norm(x1)
                    a,b,c=x1/r
                    utheta=np.array([-b,a,0.])
                    ur=np.array([a,b,0.])
                    '''
                    v=norm(v1)
                    u=np.vdot(v1,utheta)*utheta/v
                    '''
                    #--------------------------
                    v1 = np.add(v1,deltav*utheta)
                    y0 = np.hstack((x1,v1))
                if i>tendi:
                    ti = TFINAL3/nc
                    t = np.linspace(i,i+ti,nt)
                    
                    
                    #Position Terre
                    (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
                    x1 = sol1[:,0:3].T
                    v1 = sol1[:,3:6].T
                    if OrbI=='Oui':
                        plt.plot(x1[0,:],x1[1,:],'b')
                    if i<(tendi+Tterre):
                        X1 = x1.tolist()
                        TERRE[0].extend(X1[0])
                        TERRE[1].extend(X1[1])
                        TERRE[2].extend(X1[2])
                        
                    x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
                    v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
                    y1 = np.hstack((x2,v2))
        
        
                    sol2 = odeint(f,y0,t,args=(mu,))
                    x0 = sol2[:,0:3].T
                    v0 = sol2[:,3:6].T
                    if OrbI=='Oui':
                        plt.plot(x0[0,:],x0[1,:],'g')
        
                    if i<(tendi+Tasteroide):
                        X0 = x0.tolist()
                        ASTEROIDE[0].extend(X0[0])
                        ASTEROIDE[1].extend(X0[1])
                        ASTEROIDE[2].extend(X0[2])
        
                    x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                    v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                    y0 = np.hstack((x1,v1))
        
        
                    if i<(tendi+Tasteroide):
                        X0 = x0.tolist()
                        ASTEROIDE[0].extend(X0[0])
                        ASTEROIDE[1].extend(X0[1])
                        ASTEROIDE[2].extend(X0[2])
                    if etude5=='Oui':
                        periode=5
                        if i>(tendi+(periode-1)*Tasteroide) and i<(tendi+periode*Tasteroide):
                            X0 = x0.tolist()
                            ASTEROIDE3[0].extend(X0[0])
                            ASTEROIDE3[1].extend(X0[1])
                            ASTEROIDE3[2].extend(X0[2])
                        
                    x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                    v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                    y0 = np.hstack((x1,v1))
            plt.axis('equal')     
            plt.grid(True)
            plt.show()
            #Eviter les doublons
            TERRE2 = [[],[],[]] 
            for i in TERRE[0] : 
                if i not in TERRE2[0]: 
                    TERRE2[0].append(i) 
                    
            for i in TERRE[1] : 
                if i not in TERRE2[1]: 
                    TERRE2[1].append(i) 
            for i in TERRE[2] : 
                if i not in TERRE2[2]: 
                    TERRE2[2].append(i) 
            ASTEROIDE2 = [[],[],[]] 
            for i in ASTEROIDE[0] : 
                if i not in ASTEROIDE2[0]: 
                    ASTEROIDE2[0].append(i) 
            for i in ASTEROIDE[1] : 
                if i not in ASTEROIDE2[1]: 
                    ASTEROIDE2[1].append(i) 
            for i in ASTEROIDE[2] : 
                if i not in ASTEROIDE2[2]: 
                    ASTEROIDE2[2].append(i)
            if etude5=='Oui':
                ASTEROIDE4 = [[],[],[]] 
                for i in ASTEROIDE3[0] : 
                    if i not in ASTEROIDE4[0]: 
                        ASTEROIDE4[0].append(i) 
                for i in ASTEROIDE3[1] : 
                    if i not in ASTEROIDE4[1]: 
                        ASTEROIDE4[1].append(i) 
                for i in ASTEROIDE3[2] : 
                    if i not in ASTEROIDE4[2]: 
                        ASTEROIDE4[2].append(i)
            parametre=0
            
            def distmin(A):
                a=A
                L=[]
                for i in range(1,len(a),1):
                    L.append(abs(a[i]-a[i-1]))
                return(min(L))
            def distmax(A):
                a=A
                L=[]
                for i in range(1,len(a),1):
                    L.append(abs(a[i]-a[i-1]))
                return(max(L))
            print('Minimum intervalle :', distmin(TERRE2[0]),'Maximal intervalle :', distmax(TERRE2[0]))
            
            print('Nombre de points :','Terre:', len(TERRE2[0]),'Astéroide:',len(ASTEROIDE2[0]))
            if etude5=='Oui':
                print('Nombre de points Etude 5:','Nombre de périodes:',periode,'Terre:',len(TERRE2[0]),'Astéroide:',len(ASTEROIDE4[0]))
            for x in TERRE2[0]:
                for y in ASTEROIDE2[0]:
                    if abs(x-y)<10**(critere):
                        #print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                        p=TERRE2[1][TERRE2[0].index(x)]
                        c=ASTEROIDE2[1][ASTEROIDE2[0].index(y)]
                        #print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                        
                        if abs(p-c)<10**(critere):
                            '''
                            print('Test actuel :',deltav)
                            print('CRASH !!')
                            print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                            print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                            '''
                            ''' #Principe : Prendre la première valeur
                            if x<0 and P[2][indice]==0:
                                P[0][indice]=x
                                P[4][indice]=deltav
                            if p>0 and P[3][indice]==0:
                                P[1][indice]=p
                            '''
                            if x<0:
                                A[0].append(x)
                                
                            if p>0:
                                A[1].append(p)
                            #plt.plot(x,p,'o')
                            #plt.plot(y,c,'o')
                            parametre+=1
            if etude5=='Oui':                
                parametre2=0
                for x in TERRE2[0]:
                    for y in ASTEROIDE4[0]:
                        if abs(x-y)<10**(critereetude):
                            #print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                            p=TERRE2[1][TERRE2[0].index(x)]
                            c=ASTEROIDE4[1][ASTEROIDE4[0].index(y)]
                            #print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                            
                            if abs(p-c)<10**(critereetude):
                                
                                if x<0:
                                    A2[0].append(x)
                                    
                                if p>0:
                                    A2[1].append(p)
                                #plt.plot(x,p,'o')
                                #plt.plot(y,c,'o')
                                parametre2+=1
            P[0][indice]=moy(A[0])
            P[1][indice]=moy(A[1])
            P[4][indice]=deltav
            
            if etude5=='Oui':
                P2[0][indice]=moy(A2[0])
                P2[1][indice]=moy(A2[1])
                P2[4][indice]=deltav
            indice+=1
            print('Etat actuel:',indice,'Goal',precision)
            if parametre==0 and deltav!=0:
                print('Test actuel :',deltav,'Pas de crash')
                PARAMETRE[0].append(parametre)
                PARAMETRE[1].append(deltav)
                
            if etude5=='Oui':
                if parametre2==0 and deltav!=0:
                    print('Test actuel après plusieurs périodes:',deltav,'Pas de crash')
                    PARAMETRE2[0].append(parametre2)
                    PARAMETRE2[1].append(deltav)
            U[0].append(parametre)
            U[1].append(deltav)
            
         
        #plt.plot(TERRE2[0],TERRE2[1],label='Orbite de la Terre')
        #plt.plot(P[0],P[1],label='Déplacement du crash 1ere periode')
        if etude5=='Oui':
            plt.plot(P2[0],P2[1],label='Déplacement du crash derniere periode')
        
        print('1ère période',PARAMETRE[0][0],PARAMETRE[1][0])
        print(PARAMETRE)
        if etude5=='Oui':
            print('Dernière période',PARAMETRE2[0][0],PARAMETRE2[1][0])
    
    def modif(discret,dmin,dmax,ndv):
        #Etude : Programme de calcul de l'effet de modifications 
        np.seterr(divide='ignore', invalid='ignore')
        #Conditions initiales
        
        r1=0.89
        r2=2.94
        #Astéroïde
        
        r0 = max(r1,r2)*au
        a = ((r1+r2)/2)*au
        e = (abs(r1-r2))/(r1+r2)
        v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
        #print(v0)
        
        #Initialisation
        Tasteroide = np.sqrt((a**3*4*np.pi**2)/mu)
        x0 = np.array([r0,0.,0.],dtype='float64')
        v0 = np.array([0.,v0,0.],dtype='float64')
        y0 = np.hstack((x0,v0))
        
        #Terre
        r = 152*10**9
        v1 = 29.78*10**3
        Tterre = np.sqrt((r**3*4*np.pi**2)/mu)
        
        #Initialisation
        x2 = np.array([r,0.,0.],dtype='float64')
        v2 = np.array([0.,v1,0.],dtype='float64')
        y1 = np.hstack((x2,v2))
        
        #Echelle de Période
        #PERIODE=np.linspace(1,periodemax,nbp)
        PARAMETRE=[[],[]]
        
        #for periode in PERIODE:
        #print('Periode actuel :',periode)
        #Temps
        to=0
        t0=Tasteroide
        
        nb=discret
        
        impulsion=30*jour
        tiimpulsion=t0-impulsion/2
        tfinal1=np.linspace(to,tiimpulsion,nb)       # Timeline impulsion
        
        TFINAL1=tiimpulsion-to
        pas=TFINAL1/nb
        
        tendi=t0+impulsion/2
        TFINAL2=tendi-(tiimpulsion+pas)
        na=int(TFINAL2/pas)
        tfinal2=np.linspace(tiimpulsion+pas,tendi,na)
        
        tend=tendi+2*t0
        TFINAL3=tend-(tendi+pas)
        nc=int(TFINAL3/pas)
        tfinal3=np.linspace(tendi+pas,tend,nc)
        tfinal4=np.concatenate((tfinal1,tfinal2,tfinal3))
        
        #Données du protocole:
        nt=200
        critere=7
        
        

        vfinal = np.linspace(dmin,dmax,ndv)
        
        
        U=[[],[]]
        P=[[0 for i in range(ndv)] for j in range(5)]
        
        
        indice=0
        for deltav in vfinal:
            TERRE=[[],[],[]]
            ASTEROIDE = [[],[],[]]
            ASTEROIDE3 = [[],[],[]]
            A=[[],[]]
            
            for i in tfinal4:
                if i<tiimpulsion:
                    ti = pas
                    t = np.linspace(i,i+ti,nt)
                    
                    #Position Terre
                    (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
                    x1 = sol1[:,0:3].T
                    v1 = sol1[:,3:6].T
                    if OrbG=='Oui':
                        plt.plot(x1[0,:],x1[1,:],'b')
                    
                    x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
                    v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
                    y1 = np.hstack((x2,v2))
                    
                    #Position Astéroide
                    sol = odeint(f,y0,t,args=(mu,))
                    x0 = sol[:,0:3].T
                    v0 = sol[:,3:6].T
                    if OrbG=='Oui':
                        plt.plot(x0[0,:],x0[1,:],'g')
                    x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
                    v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
                    y0 = np.hstack((x1,v1))
            
                    
                if i>=tiimpulsion and i<=tendi: #échelle de temps de l'impulsion
                    ti = pas
                    t = np.linspace(i,i+ti,nt)
                    
                    #Position Terre
                    (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
                    x1 = sol1[:,0:3].T
                    v1 = sol1[:,3:6].T
                    if OrbH=='Oui':
                        plt.plot(x1[0,:],x1[1,:],'b')
                    
                    x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
                    v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
                    y1 = np.hstack((x2,v2))
                    
                    #Position Astéroide
                    sol = odeint(f,y0,t,args=(mu,))
                    x0 = sol[:,0:3].T
                    v0 = sol[:,3:6].T
                    if OrbH=='Oui':
                        plt.plot(x0[0,:],x0[1,:],'g')
                    
                    x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                    v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                    #-------------------
                    r=norm(x1)
                    a,b,c=x1/r #a=cos(theta), b=sin(theta)
                    utheta=np.array([-b,a,0.])
                    '''
                    v=norm(v1)
                    u=np.vdot(v1,utheta)*utheta/v
                    '''
                    #--------------------------
                    v1 = np.add(v1,deltav*utheta)
                    y0 = np.hstack((x1,v1))
                if i>tendi:
                    ti = pas
                    t = np.linspace(i,i+ti,nt)
                    
                    
                    #Position Terre
                    (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
                    x1 = sol1[:,0:3].T
                    v1 = sol1[:,3:6].T
                    if OrbI=='Oui':
                        plt.plot(x1[0,:],x1[1,:],'b')
                    if i<(tendi+Tterre):
                        X1 = x1.tolist()
                        TERRE[0].extend(X1[0])
                        TERRE[1].extend(X1[1])
                        TERRE[2].extend(X1[2])
                        
                    x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
                    v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
                    y1 = np.hstack((x2,v2))
        
        
                    sol2 = odeint(f,y0,t,args=(mu,))
                    x0 = sol2[:,0:3].T
                    v0 = sol2[:,3:6].T
                    if OrbI=='Oui':
                        plt.plot(x0[0,:],x0[1,:],'g')
        
                    if i<(tendi+Tasteroide):
                        X0 = x0.tolist()
                        ASTEROIDE[0].extend(X0[0])
                        ASTEROIDE[1].extend(X0[1])
                        ASTEROIDE[2].extend(X0[2])
        
                    x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                    v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                    y0 = np.hstack((x1,v1))
        
        
                    if i<(tendi+Tasteroide):
                        X0 = x0.tolist()
                        ASTEROIDE[0].extend(X0[0])
                        ASTEROIDE[1].extend(X0[1])
                        ASTEROIDE[2].extend(X0[2])
                    if etude5=='Oui':
                        periode=5
                        if i>(tendi+(periode-1)*Tasteroide) and i<(tendi+periode*Tasteroide):
                            X0 = x0.tolist()
                            ASTEROIDE3[0].extend(X0[0])
                            ASTEROIDE3[1].extend(X0[1])
                            ASTEROIDE3[2].extend(X0[2])
                        
                    x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
                    v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
                    y0 = np.hstack((x1,v1))
            #Eviter les doublons
            TERRE2 = [[],[],[]] 
            for i in TERRE[0] : 
                if i not in TERRE2[0]: 
                    TERRE2[0].append(i) 
                    
            for i in TERRE[1] : 
                if i not in TERRE2[1]: 
                    TERRE2[1].append(i) 
            for i in TERRE[2] : 
                if i not in TERRE2[2]: 
                    TERRE2[2].append(i) 
            ASTEROIDE2 = [[],[],[]] 
            for i in ASTEROIDE[0] : 
                if i not in ASTEROIDE2[0]: 
                    ASTEROIDE2[0].append(i) 
            for i in ASTEROIDE[1] : 
                if i not in ASTEROIDE2[1]: 
                    ASTEROIDE2[1].append(i) 
            for i in ASTEROIDE[2] : 
                if i not in ASTEROIDE2[2]: 
                    ASTEROIDE2[2].append(i)
            if (0)==0:
                #Brouillon
                parametre=0
                
            
                if indice==0:
                    print('Discretisation:',discret,'Nombre de points :','Terre:', len(TERRE2[0]),'Astéroide:',len(ASTEROIDE2[0]))
                
                #Tri rapide
                def trirap(a):
                    if len(a)==0:
                        return(a)
                    a1,a2=[],[]
                    for x in a[1:]:
                        if x<=a[0]:
                            a1.append(x)
                        else:
                            a2.append(x)
                    return(trirap(a1)+[a[0]]+trirap(a2))
                #Calcul de la moyenne de la distance entre deux points dans une liste
                def distmoy(A):
                    #a=trirap(A)
                    a=A
                    L=[]
                    for i in range(1,len(a),1):
                        L.append(abs(a[i]-a[i-1]))
                    return('Moyenne :',moy(L),'Minimum :',min(L))
                '''
                if indice==0:
                    print(distmoy(TERRE2[0]))
                '''
                for x in TERRE2[0]:
                    for y in ASTEROIDE2[0]:
                        if abs(x-y)<10**(critere):
                            #print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                            p=TERRE2[1][TERRE2[0].index(x)]
                            c=ASTEROIDE2[1][ASTEROIDE2[0].index(y)]
                            #print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                            
                            if abs(p-c)<10**(critere):
                                
                                '''
                                print('Test actuel :',deltav)
                                print('CRASH !!')
                                print('Coord. X Terre:',x,'Coord. X Asteroide :',y)
                                print('Coord Y Terre :',p,'Coord. Y Asteroide :',c)
                                '''
                                
                                
                                '''
                                #Principe : Prendre la première valeur
                                if x<0 and P[2][indice]==0:
                                    P[0][indice]=x
                                    P[4][indice]=deltav
                                if p>0 and P[3][indice]==0:
                                    P[1][indice]=p
                                '''
    
                                if x<0:
                                    A[0].append(x)
                                    
                                if p>0:
                                    A[1].append(p)
                                #plt.plot(x,p,'o')
                                #plt.plot(y,c,'o')
                                parametre+=1
                
                P[0][indice]=moy(A[0])
                P[1][indice]=moy(A[1])
                P[4][indice]=deltav
                
                
                indice+=1
                print('Progression :',indice/ndv*100)
                #print('Etat actuel:',indice,'Goal',precision)
                PARAMETRE[0].append(parametre)
                PARAMETRE[1].append(deltav)
                
                    
                
                U[0].append(parametre)
                print(parametre)
                U[1].append(deltav)
                
            
            #plt.plot(TERRE2[0],TERRE2[1],label='Orbite de la Terre')
            #plt.plot(P[0],P[1],label='Déplacement du crash 1ere periode')
            
            
            #print('1ère période',PARAMETRE[0][0],PARAMETRE[1][0])
        return(PARAMETRE)
            
        '''
            print('Discretisation:',discret,'Nombre de points :','Terre:', len(TERRE2[0]),'Astéroide:',len(ASTEROIDE2[0]))
            INTERSECTION=[[],[]]
            #Méthode cKDTree
            from scipy.spatial import cKDTree
            from scipy import interpolate
            
            
            
            def upsample_coords(coord_list):
            # s is smoothness, set to zero
            # k is degree of the spline. setting to 1 for linear spline
                tck, u = interpolate.splprep(coord_list, k=5, s=0.0)
                upsampled_coords = interpolate.splev(np.linspace(0, 1, 100), tck)
                return upsampled_coords
            
            # target line
            
            
            x_targ = TERRE2[0]
            y_targ = TERRE2[1]
            z_targ = [0 for i in range(len(TERRE2[0]))]
            targ_upsampled = upsample_coords([x_targ, y_targ, z_targ])
            targ_coords = np.column_stack(targ_upsampled)
            # KD-tree for nearest neighbor search
            targ_kdtree = cKDTree(targ_coords)
            
            # line two
            x2 = ASTEROIDE2[0]
            y2 = ASTEROIDE2[1]
            z2 = [0 for i in range(len(ASTEROIDE2[0]))]
            l2_upsampled = upsample_coords([x2, y2, z2])
            l2_coords = np.column_stack(l2_upsampled)
            
            # find intersections
            for i in range(len(l2_coords)):
                if i == 0:  # skip first, there is no previous point
                    continue
                distance, close_index = targ_kdtree.query(l2_coords[i], distance_upper_bound=10**(9))
            # strangely, points infinitely far away are somehow within the upper bound
                if np.isinf(distance):
                    continue
            # plot ground truth that was activated
                _x, _y, _z = targ_kdtree.data[close_index]
                INTERSECTION[0].append([_x,_y,_z])
                
                _x2, _y2, _z2 = l2_coords[i]
                INTERSECTION[1].append([_x2,_y2,_z2])
            U[0].append(len(INTERSECTION[0]))
            U[1].append(deltav)
        return U
            '''
        
    if etude7=='Oui':
        #Recherche du point pivot par dichotomie en considérant la fonction modif
        ETUDE7=[[],[],[]]
        etude7liste=np.linspace(100,300,10)
        for discret in etude7liste:
            ETUDE7[0].append(int(discret))
            a=0
            b=10000
            while(b-a)>10**(-1):
                PARAMETRE=modif(int(discret),a,b,5)
                compteur=0
                for p in PARAMETRE[0]:
                    if compteur==0 and p==0:
                        a=PARAMETRE[1][PARAMETRE[0].index(p)-1]
                        b=PARAMETRE[1][PARAMETRE[0].index(p)]
                        compteur+=1
                        print('discret:','a=',a,'b=',b)
            ETUDE7[1].append(a)
            ETUDE7[2].append(b)
        import statistics
        v1=statistics.pstdev(ETUDE7[1], mu=None)
        v2=statistics.pstdev(ETUDE7[2], mu=None)
        m1=moy(ETUDE7[1])
        m2=moy(ETUDE7[2])
        print('Moy 1',m1,'Moy2',m2,'Var1',v1,'Var2',v2)
            
#-----------------------------------------------   
if protocole6=='Activé':
    #Protocole 6 = Ouverture à la 3D
    
    
 
    np.seterr(divide='ignore', invalid='ignore')
    #Conditions initiales
    
    #Astéroïde
    r1=0.89
    r2=2.94
    inclinaison=10 #Inclinaison par rapport au plan de l'éclipitique en degrés 
    
    r0 = max(r1,r2)*au
    a = ((r1+r2)/2)*au
    e = (abs(r1-r2))/(r1+r2)
    v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
    inclinaison=inclinaison*np.pi/180
    
    #Initialisation
    Tasteroide = np.sqrt((a**3*4*np.pi**2)/mu)
    x0 = np.array([r0*np.cos(inclinaison),0.,r0*np.sin(inclinaison)],dtype='float64')
    v0 = np.array([0.,v0,0.],dtype='float64')
    y0 = np.hstack((x0,v0))
    
    #Terre
    r = 152*10**9
    v1 = 29.78*10**3
    Tterre = np.sqrt((r**3*4*np.pi**2)/mu)

    #Initialisation
    x2 = np.array([r,0.,0.],dtype='float64')
    v2 = np.array([0.,v1,0.],dtype='float64')
    y1 = np.hstack((x2,v2))
    
    #Temps
    to=0
    t0=Tasteroide
    
    nb = 1000
    impulsion=30*jour
    tiimpulsion=t0-impulsion/2
    
    tfinal1=np.linspace(to,tiimpulsion,nb)       # Timeline impulsion
    TFINAL1=tiimpulsion-to
    pas=TFINAL1/nb
    
    tendi=t0+impulsion/2
    TFINAL2=tendi-(tiimpulsion+pas)
    na=int(TFINAL2/(0.1*jour))
    
    tfinal2=np.linspace(tiimpulsion+pas,tendi,na)
    tend=tendi+2*t0
    TFINAL3=tend-(tendi+pas)
    nc = 40000
    tfinal3=np.linspace(tendi+pas,tend,nc)
    tfinal4=np.concatenate((tfinal1,tfinal2,tfinal3))
    n = len(tfinal4)
    
    #Données du protocole:
    nt=2
    
    #Plots
    #fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot([0.],[0.],'o',label='Centre attracteur')
    #Multiplication d'un scalaire sur une liste
    def mult(a, A):
        L = []
        for i in range(len(A)):
            L.append(A[i] * a)
        return L
    #Soustraction de deux listes élément par élément
    def sous(A,B):
        L = []
        for i in range(len(A)):
            L.append(A[i] - B[i])
        return L
    
    
    def g(y1,t,y2,mu,mp):
        x1 = y1[0:3].copy()
        x2 = y2[0:3].copy()
        v1 = y1[3:6].copy()
        r = np.sqrt(x1[0]**2+x1[1]**2+x1[2]**2)
        p = np.sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2+(x1[2]-x2[2])**2)
        dxdt = v1
        dvdt =  (mult(-mu,x1))/(r**3) + (mult(-mp,sous(x1,x2)))/(p**3) 
        dy = np.hstack((dxdt,dvdt))
        return dy
    
    mplanete = 6*10**24
    mp = 6.67*10**(-11)*mplanete                  # Gravitational parameter planet
    TERRE=[[],[],[]]
    ASTEROIDE = [[],[],[]]
    for i in tfinal4:
        if i<=t0:
            ti = TFINAL1/nb
            t = np.linspace(i,i+ti,nt)
            
            ''' #Sans prendre en compte l'attraction terrestre
            
            sol = odeint(f,y0,t,args=(mu,))
            x0 = sol[:,0:3].T
            v0 = sol[:,3:6].T
            #orb1=plt.plot(x0[0,:],x0[1,:],'g')
            x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
            v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
            y0 = np.hstack((x1,v1))
            '''
            #Position Terre
            (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
            x1 = sol1[:,0:3].T
            v1 = sol1[:,3:6].T
            #plt.plot(x1[0,:],x1[1,:],'b')
            
                
            x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
            v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
            y1 = np.hstack((x2,v2))


            (sol2,d) = odeint(g,y0,t,args=(y1,mu,mp,),full_output=1)
            x0 = sol2[:,0:3].T
            v0 = sol2[:,3:6].T
            #plt.plot(x0[0,:],x0[1,:],'g')

            

            x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
            v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
            y0 = np.hstack((x1,v1))
            
        if i>t0 and i<=tendi: #échelle de temps de l'impulsion
            ti = TFINAL2/na
            t = np.linspace(i,i+ti,nt)
            
            ''' #Sans prendre en compte l'attraction terrestre
            sol = odeint(f,y0,t,args=(mu,))
            x0 = sol[:,0:3].T
            v0 = sol[:,3:6].T
            #orb1=plt.plot(x0[0,:],x0[1,:],'r')
            x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
            v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
            v1 = np.add(v1,deltav*u)                             # IMPULSION
            y0 = np.hstack((x1,v1))
           '''
           
           #Position Terre
            (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
            x1 = sol1[:,0:3].T
            v1 = sol1[:,3:6].T
            #plt.plot(x1[0,:],x1[1,:],'b')
            
                
            x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
            v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
            y1 = np.hstack((x2,v2))


            (sol2,d) = odeint(g,y0,t,args=(y1,mu,mp,),full_output=1)
            x0 = sol2[:,0:3].T
            v0 = sol2[:,3:6].T
            #plt.plot(x0[0,:],x0[1,:],'g')

            

            x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
            v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
            #v1 = np.add(v1,deltav*u)
            y0 = np.hstack((x1,v1))
        if i>tendi:
            ti = TFINAL3/nc
            t = np.linspace(i,i+ti,nt)
            
            ''' #Sans prendre en compte l'attraction terrestre
            sol = odeint(f,y0,t,args=(mu,))
            x0 = sol[:,0:3].T
            v0 = sol[:,3:6].T
            #orb1=plt.plot(x0[0,:],x0[1,:],'g',label='Après impulsion discrétisée')
            x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]])
            v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]])
            y0 = np.hstack((x1,v1))   
            '''
            #Position Terre
            (sol1,d) = odeint(f,y1,t,args=(mu,),full_output=1)
            x1 = sol1[:,0:3].T
            v1 = sol1[:,3:6].T
            #plt.plot(x1[0,:],x1[1,:],'b')
            if i<(tendi+Tterre):
                X1 = x1.tolist()
                TERRE[0].extend(X1[0])
                TERRE[1].extend(X1[1])
                TERRE[2].extend(X1[2])
                
            x2 = np.array([x1[0][nt-1],x1[1][nt-1],x1[2][nt-1]])
            v2 = np.array([v1[0][nt-1],v1[1][nt-1],v1[2][nt-1]])
            y1 = np.hstack((x2,v2))


            (sol2,d) = odeint(g,y0,t,args=(y1,mu,mp,),full_output=1)
            x0 = sol2[:,0:3].T
            v0 = sol2[:,3:6].T
            #plt.plot(x0[0,:],x0[1,:],'g')

            if i<(tendi+Tasteroide):
                X0 = x0.tolist()
                ASTEROIDE[0].extend(X0[0])
                ASTEROIDE[1].extend(X0[1])
                ASTEROIDE[2].extend(X0[2])

            x1 = np.array([x0[0][nt-1],x0[1][nt-1],x0[2][nt-1]]) # Récursivité
            v1 = np.array([v0[0][nt-1],v0[1][nt-1],v0[2][nt-1]]) # Récursivité
            y0 = np.hstack((x1,v1))
    #Eviter les doublons
    TERRE2 = [[],[],[]] 
    for i in TERRE[0] : 
        if i not in TERRE2[0]: 
            TERRE2[0].append(i) 
    for i in TERRE[1] : 
        if i not in TERRE2[1]: 
            TERRE2[1].append(i) 
    for i in TERRE[2] : 
        if i not in TERRE2[2]: 
            TERRE2[2].append(i) 
    ASTEROIDE2 = [[],[],[]] 
    for i in ASTEROIDE[0] : 
        if i not in ASTEROIDE2[0]: 
            ASTEROIDE2[0].append(i) 
    for i in ASTEROIDE[1] : 
        if i not in ASTEROIDE2[1]: 
            ASTEROIDE2[1].append(i) 
    for i in ASTEROIDE[2] : 
        if i not in ASTEROIDE2[2]: 
            ASTEROIDE2[2].append(i)
    
    ax.plot3D(TERRE2[0], TERRE2[1], TERRE2[2], 'red',label="Orbite Terre")
    ax.plot3D(ASTEROIDE2[0], ASTEROIDE2[1], ASTEROIDE2[2], 'blue',label="Astéroide")
    #plt.axis('equal')     # ¡ to see circular orbits as a circle !
    #plt.grid(True)
    plt.title('Orbites and co. Application à 2019 PDC')
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.xlabel("x")
    plt.ylabel("y")
    #plt.zlabel("z")
    
    
 
#----------------------------------------------- 
if protocole8=='Activé':
    #Protocole 8 = Simulation de plusieurs astéroïdes et du Système solaire
    if asteroidonly=="Oui":
        #Astéroïde
        nasteroide=10
        for i in range(nasteroide):
            r1=6*np.random.random_sample(1)
            r2=6*np.random.random_sample(1)
            
            r0 = max(r1,r2)*au
            a = ((r1+r2)/2)*au
            e = (abs(r1-r2))/(r1+r2)
            v0 = np.sqrt((mu*(1-e)/(a*(1+e))))
            
            x0 = np.array([r0,0.,0.],dtype='float64')
            v0 = np.array([0.,v0,0.],dtype='float64')
            y0 = np.hstack((x0,v0))
            
            to=0
            tend=10*365*jour
            nt=1500
            tfinal=np.linspace(to,tend,nt)
            sol = odeint(f,y0,tfinal,args=(mu,))
            x0 = sol[:,0:3].T
            orb2=plt.plot(x0[0,:],x0[1,:],'c')
    if solarsystem=='Oui':
        #Cf. Problème à N corps - Euler.py
        print('Cf. Programme extérieur')
#-----------------------------------------------   

if protocoleTerre=='Activé':
    #Protocole Terre : Orbite de la Terre
    #Paramètres de la Terre
    r0 = 152*10**9
    v0 = 29.78*10**3
    x0 = np.array([r0,0.,0.],dtype='float64')
    v0 = np.array([0.,v0,0.],dtype='float64')
    y0 = np.hstack((x0,v0))
    nt = 5000
    
    to=0
    tend=1.5*365*jour
    #Orbite de base
    tfinal=np.linspace(to,tend,nt)
    sol = odeint(f,y0,tfinal,args=(mu,))
    x0 = sol[:,0:3].T
    v0 = sol[:,3:6].T
    orb2=plt.plot(x0[0,:],x0[1,:],'c',label='Orbite terrestre')

#-----------------------------------------------
if plot=='Oui':
    #Plots
    plt.plot([0.],[0.],'o',label='Centre attracteur')
    plt.axis('equal')     # ¡ to see circular orbits as a circle !
    plt.grid(True)
    plt.title('Collisions')
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
#End Orbit-propagation-with-odeint (Modified)
