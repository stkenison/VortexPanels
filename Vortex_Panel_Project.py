import numpy as np; import os;import matplotlib.pyplot as plt; import json
os.system('cls'); plt.rcParams['font.family'] = 'Times New Roman'

#get inputs from .json file
with open('input.json') as f: data = json.load(f)

#initialize variables
geometry = np.loadtxt(data['geometry'],dtype = float) #Load airfoil geometry from file
alpha = np.array(data['alpha[deg]'])/180*np.pi #define desired angles of attack
v_inf = data['freestream_velocity']
n = geometry.shape[0] #count number of points in airfoil geometry file
n_alpha = len(alpha) #count number of unique values in angle of attack list

#initialize variables to describe geometry
CP = 0.5*(geometry[1:]+geometry[:-1])#define control points
dx = geometry[1:,0]-geometry[:-1,0] #define dx values
dy = geometry[1:,1]-geometry[:-1,1] #define dy values
l = np.sqrt(dx**2+dy**2) #define lengths of panels

#solve for A array
A = np.zeros((n,n)) #initialize A array
for i in range(n-1):
    for j in range(n-1):
        xi = (dx[j]*(CP[i,0]-geometry[j,0])+dy[j]*(CP[i,1]-geometry[j,1]))/l[j] #define xi term
        eta = (-dy[j]*(CP[i,0]-geometry[j,0])+dx[j]*(CP[i,1]-geometry[j,1]))/l[j] #define eta term
        
        phi = np.arctan2(eta*l[j],eta**2+xi**2-xi*l[j]) #define phi term
        psi = 0.5*np.log((xi**2+eta**2)/((xi-l[j])**2+eta**2)) #define psi term

        P1 = np.array([[dx[j],-dy[j]],[dy[j],dx[j]]]) #calculate first part of P array
        P2 = np.array([[(l[j]-xi)*phi+eta*psi,xi*phi-eta*psi],
            [eta*phi-(l[j]-xi)*psi-l[j],-eta*phi-xi*psi+l[j]]]) #calculate second part of P array
        P = (1/(2*np.pi*l[j]**2))*np.matmul(P1,P2) #create P array by matrix multiplying parts

        A[i,j] += dx[i]*P[1,0]/l[i]-dy[i]*P[0,0]/l[i] #calculat current cell of A array
        A[i,j+1] += dx[i]*P[1,1]/l[i]-dy[i]*P[0,1]/l[i] #start calculating next cell of A array

A[n-1,0] = 1.0 #define edge of A array based on kutta condition
A[n-1,n-1] = 1.0 #define corner of A array based on kutta condition

#use B array to solve for gamma for all angles of attack
B = np.zeros((n,n_alpha)) #initialize B column(s)
gamma = np.zeros((n,n_alpha))  #initialize gamma column(s)
for i in range(n_alpha):
    B[:-1,i] = v_inf*((geometry[1:,1]-geometry[:-1,1])*np.cos(alpha[i])-
        (geometry[1:,0]-geometry[:-1,0])*np.sin(alpha[i]))/[l[:]] #create B array for current angles of attack
    gamma[:,i] = np.linalg.solve(A,B[:,i]) #solve for gamma asoociated with current angle of attack

#find coeffecient of lift
C_L = np.zeros(n_alpha) #initialize lift coefficient array
for i in range(n_alpha):
    C_L[i] = np.sum(l[:]*(gamma[:-1,i]+gamma[1:,i])/v_inf) #find value of current iteration of lift coefficient array

#find leading edge moment coefficient array
C_mLE = np.zeros(n_alpha) #initialize leading edge moment coefficient array
for i in range(n_alpha):
    C_mLE[i] = -1/3*np.sum(l[:]/v_inf*((2*geometry[:-1,0]*gamma[:-1,i]+geometry[1:,0]*gamma[:-1,i]+geometry[:-1,0]*gamma[1:,i]+2*geometry[1:,0]*gamma[1:,i])*np.cos(alpha[i])
    +(2*geometry[:-1,1]*gamma[:-1,i]+geometry[1:,1]*gamma[:-1,i]+geometry[:-1,1]*gamma[1:,i]+2*geometry[1:,1]*gamma[1:,i])*np.sin(alpha[i]))) #find value of current iteration of lift coefficient array

# Find quarter chord moment coefficient array
C_mQC = np.zeros(n_alpha)  #initialize quarter chord moment coefficient array
QC_x = geometry[:,0]-0.25 #Compute the moment arm relative to the quarter chord
for i in range(n_alpha):
    C_mQC[i] = -1/3 * np.sum(
        l[:]/v_inf*(
            (2*QC_x[:-1]*gamma[:-1, i]+QC_x[:-1]*gamma[1:,i]+ 
            QC_x[1:]*gamma[:-1,i]+2*QC_x[1:]*gamma[1:,i])*np.cos(alpha[i])+
            (2*geometry[:-1,1]*gamma[:-1,i]+geometry[1:,1]*gamma[:-1,i]+ 
            geometry[:-1,1]*gamma[1:,i]+2*geometry[1:,1]*gamma[1:,i])*np.sin(alpha[i])
        )
    )  # find value of current iteration of quarter chord moment coefficient array

#print out calculated coefficients for verification
print(f"C_L: {C_L}")
print(f"C_mLE: {C_mLE}")
print(f"C_mQC: {C_mQC}")


#modify to read in inputs from json file
#modify code to loop through various angles of attack

plt.figure(1)
plt.plot(geometry[:,0],geometry[:,1])
plt.axis("equal")
plt.title('Airfoil Geometry')

plt.figure(2)
plt.plot(alpha*180/np.pi,C_L)
plt.title('Lift Coefficient ($C_{L}$)')
plt.xlabel("Angle ($^{o}$)");plt.ylabel("$C_{L}$",rotation='horizontal')
#plt.show()
#plt.savefig('Figure{}'.format(3+k*3));plt.pause(0.00001)
