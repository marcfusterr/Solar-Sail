#Creator Marc Fuster, Daniel Goncalveç and Nicetu Tibau

#This code was created for a 48h competition. There was not any time. Therefore
# the code is not clean nor commented.

from numpy.linalg import norm as norm
import numpy as np
import matplotlib.pyplot as plt

#simple arctan function
def arctan(y,x):
    if x==0:
        a=np.pi/2*float(y)/np.abs(y)
    else:
        a=np.arctan(float(y)/x)+np.pi*(np.abs(x)-x)/(2*np.abs(x))
    return a

def accelerationEarth(rt,rm):
    a=-1*G*(Ms*rt/(norm(rt)**3)+MM*(rt-rm)/(norm(rt-rm)**3))
    return a

def accelerationMars(rt,rm):
    a=-G*(Ms*rm/(norm(rm)**3)+Mt*(rm-rt)/(norm(rm-rt)**3))
    return a

def accelerationShip(rt,rm,r,Psi_s,Psi_m,Theta,A):
    a_g=-G*(Ms*r/(norm(rm)**3)+MM*(rm-r)/(norm(rm-r)**3)+Mt*(rt-r)/(norm(rt-r))**3)
    ax_l=0.#acceleration due to radiation pressure in x direction
    ay_l=0.#acceleration due to radiation pressure in y direction
    CE=np.dot(r-rm/norm(r-rm),-r/norm(r))  
    
    #far from a planet
    if -np.pi/2 <= -Psi_s+Theta <= np.pi/2:
        ax_l=ax_l+A/m/(np.pi*c*2)*Ls/(norm(r))**2*(np.sin(Psi_s-Theta))**3*((np.cos(Psi_s-Theta))*r[0]-(np.sin(Psi_s-Theta))*r[1])/norm(r)#Sun is very far
        ay_l=ay_l+A/m/(np.pi*c*2)*Ls/(norm(r))**2*(np.sin(Psi_s-Theta))**3*((np.cos(Psi_s-Theta))*r[1]+(np.sin(Psi_s-Theta))*r[0])/norm(r)#Sun is very far

    #close to a planet
    if -np.pi/2 <= -Psi_m+Theta <= np.pi/2:
        
        CE=np.dot((r-rm),r)/norm(r)/norm(r-rm)
        ax_l=ax_l+A/m/(12*np.pi*c*norm(rm)**2)*albedo*Ls*(1.-(1.-norm(Rm)/norm(r))**(3/2.))*(1.+CE)/2*((np.cos(Psi_m-Theta))*r[0]-(np.sin(Psi_m-Theta))*r[1])/norm(r)
        ay_l=ay_l+A/m/(12*np.pi*c*norm(rm)**2)*albedo*Ls*(1.-(1.-norm(Rm)/norm(r))**(3/2.))*(1.+CE)/2*((np.cos(Psi_m-Theta))*r[1]+(np.sin(Psi_m-Theta))*r[0])/norm(r)      

    a=np.array([0.,0.])
    a[0]=a_g[0]+ax_l
    a[1]=a_g[1]+ay_l
    return a 

#Do the computations
def kernel(numberiterations,rti,rmi,ri,vti,vmi,vi,A,a):

    #to store
    rearth=[]
    rmars=[]
    vearth=[]
    vmars=[]

    rship=[]
    vship=[]
    
    #initial condition
    rt=np.array(rti)
    rm=np.array(rmi)
    r=np.array(ri)  
    vt=np.array(vti)
    vm=np.array(vmi)
    v=vi
    
    breakconditionearth=[]
    breakconditionmars=[]    
    
    h=float(10)
    
    #plt.figure()
    
    #just initial values they will be updatet
    Psi_s=0.
    Psi_m=np.pi                    
    Theta=0.2  
    
  
    #iteration loops
    for i in range(numberiterations):
                   
        k1velearth=accelerationEarth(rt,rm)
        k1velmars=accelerationMars(rt,rm)
        k1velship=accelerationShip(rt,rm,r,Psi_s,Psi_m,Theta,A)

        k2velearth=accelerationEarth(rt+h*k1velearth/2,rm+h*k1velmars/2)
        k2velmars=accelerationMars(rt+h*k1velearth/2,rm+h*k1velmars/2)
        k2velship=accelerationShip(rt+h*k1velearth/2,rm+h*k1velmars/2,r+h*k1velship/2,Psi_s,Psi_m,Theta,A)    
                
        k3velearth=accelerationEarth(rt+h*k2velearth/2,rm+h*k2velmars/2)
        k3velmars=accelerationMars(rt+h*k2velearth/2,rm+h*k2velmars/2) 
        k3velship=accelerationShip(rt+h*k2velearth/2,rm+h*k2velmars/2,r+h*k2velship/2,Psi_s,Psi_m,Theta,A) 
    
        k4velearth=accelerationEarth(rt+h*k3velearth,rm+h*k3velmars)
        k4velmars=accelerationMars(rt+h*k3velearth,rm+h*k3velmars)
        k4velship=accelerationShip(rt+h*k3velearth,rm+h*k3velmars,r+h*k3velship,Psi_s,Psi_m,Theta,A)    
    
        k1posearth=vt
        k1posmars=vm
        k1posship=v
            
        k2posearth=k1velearth*h/2+vt
        k2posmars=k1velmars*h/2+vm
        k2posship=k1velship*h/2+v    
   
        k3posearth=vt+k2velearth*h/2
        k3posmars=vm+k2velmars*h/2
        k3posship=v+k2velship*h/2    
    
        k4posearth=vt+k3velearth*h
        k4posmars=vm+k3velmars*h
        k4posship=v+k3velship*h
                
        vt=vt+h*(k1velearth+2*k2velearth+2*k3velearth+k4velearth)/6 
        vm=vm+h*(k1velmars+2*k2velmars+2*k3velmars+k4velmars)/6         
        v=v+h*(k1velship+2*k2velship+2*k3velship+k4velship)/6
    
        rt=rt+h*(k1posearth+2*k2posearth+2*k3posearth+k4posearth)/6
        rm=rm+h*(k1posmars+2*k2posmars+2*k3posmars+k4posmars)/6
        r=r+h*(k1posship+2*k2posship+2*k3posship+k4posship)/6
        
        distancetoearth=norm(r-rt)
        distancetomars=norm(r-rm)
        breakconditionearth.append(distancetoearth)  
        breakconditionmars.append(distancetomars)   
            
        rearth.append(rt)
        rmars.append(rm)
        rship.append(r)
        vearth.append(vt)    
        vmars.append(vm)
        vship.append(v)
        
        Psi_s=float(arctan(r[1],r[0]))
        Psi_m=float(arctan((r[1]-rm[1]),(r[0]-rm[0])))
        
        vsolution=12.
        rsolution=12.

        
        if((rm[0]-r[0])**2+(rm[1]-r[1])**2)**0.5<a:
            Theta=Psi_m-Psi_s
            if i%300 ==0:
                    plt.scatter(rt[0],rt[1],color='c')
                    plt.scatter(rm[0],rm[1],color='red')    
                    plt.scatter(r[0],r[1],color='green')
                    plt.pause(0.005)
            #print('frenamos')
            
               
        else:    
            Theta=Psi_s+np.pi/2
            if i%300 ==0:
                plt.scatter(rt[0],rt[1],color='c')
                plt.scatter(rm[0],rm[1],color='red')    
                plt.scatter(r[0],r[1],color='black')
                plt.pause(0.005)
        
#        if auxtheta-Theta==0 :#or auxtheta-Theta>0:
#            plt.scatter(r[0],r[1],color='red')
#            print('canvia a')
#            print(r)
                    
        
        if distancetoearth<Rt:
            print('Entered Earth')
            break
        if distancetomars<10*Rm:
            print('It entered Mars')
            print('the index is')
            print(i)
            print('the velocity is:')
            print(vm-v)
            vsolution=vm-v
            rsolution=distancetomars-Rm
            
            break
   
        if distancetoearth<specialstepdistanceearth or distancetomars<0.1*specialstepdistancemars: 
            h=float(10)
            #print('x')
        
        elif distancetomars<specialstepdistancemars:
            h=float(100)
            
        else:
            h=float(1000)        
        


       # if r[1]>0.2e10+rm[1] and rm[0]+0.2e10>r[0]>rm[0]-0.2e10:
        #    print('se paso')
         #   break
            
    minvalue=min(breakconditionmars)
    print(minvalue)
    indexminvalue=breakconditionmars.index(min(breakconditionmars))
    print(indexminvalue)
    velocityatmin=vship[indexminvalue]
    marsvelocity=vmars[indexminvalue]  
    vsolution=norm(marsvelocity-velocityatmin)
    rsolution=minvalue

    rearth=np.matrix(rearth)    
    rship=np.matrix(rship)    
    rmars=np.matrix(rmars)    
      

    return minvalue,velocityatmin,marsvelocity,rmars,indexminvalue,rship,vsolution,rsolution


#beggining of the script
plt.close('all')

G=float(6.67428e-11)
c=float(299792458)
Au=float(149597870700)

#Planets conditions
Ms=float(1.989e30)
Mt=float(5.9723e24)
MM=float(0.64171e24)
m=2000

leo=float(2e7)
Rt=float(6.371e6)
Rm=float(3.390e6)

albedo=0.25
Ls=float(3.846e26)

#numberiterations=35000
numberiterations=50000
divisions=1
#lpha=np.linspace(2.0995,2.0995,divisions)
alpha=2.69

vesc=(2*G*Mt/(Rt+leo))**0.5


vti=np.array([0,float(29780)])
vmi=np.array([0,float(26500)])
vi=np.array([vesc*np.sin(alpha),vti[1]+vesc*np.cos(alpha)])
rti=np.array([float(152.098e9),0])
rmi=np.array([float(206.655e9),0])
ri=np.array([rti[0]+Rt+leo,0])

specialstepdistanceearth=norm(rti-rmi)/500
specialstepdistancemars=norm(rti-rmi)/8

Avec=np.array([1800,1800,1900])/7e-03
acond=[norm(rti-rmi)*0.6415,norm(rti-rmi)/8,norm(rti-rmi)/16]

minimumdistance=np.zeros([3,3])
velocitiatminimum=np.zeros([3,3])
marsvelatmin=np.zeros([3,3])
solution=np.zeros([3,3])
#

print('la iteracio es')
A=Avec[0]
a=acond[0]
print(a,A)
plt.figure()
[auxdist,auxveloc,auxmarsvelocity,rmars,indexminvalue,rship,vsolution,rsolution]=kernel(numberiterations,rti,rmi,ri,vti,vmi,vi,A,a)
minimumdistance[0][0]=auxdist
    



velocitiatminimum=np.array(velocitiatminimum)
minimumdistance=np.array(minimumdistance)
marsvelatmin=np.array(marsvelatmin)
   
#plt.title(r'Orbit with $\theta=2.69$',fontsize=23)
#plt.xlabel(r'X (m)',fontsize=20)
#plt.ylabel(r'Y (m)',fontsize=20)   
#plt.scatter(0,0,marker="o")
#plt.xlim([-1e11,3e11])
#plt.ylim([-0.5e11,2e11])

















    
    
    