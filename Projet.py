import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as sp 

def Mat(a,b,K,f,N):
	h = (b-a)/(N-1)
	t = np.linspace(a,b,N)
	x = np.linspace(a,b,N) 
	A = np.ones((N,N))
	I = np.identity(N)
	F = np.zeros(N)
	for i in range (N):
		F[i] = f(x[i])
		for j in range (N):
			A[i,j] = 2*K(x[i],t[j])

	A[:,0]=A[:,0]/2
	A[:,-1]=A[:,-1]/2
	F=np.transpose(F)
	M = (I - (h/2)*A)

	return t,F,M


#3.2

def K_test(x,t):
	return(2)

def F_test(x):
	return (np.cos((np.pi*x)/2) -(8/np.pi))

def U_test(x):
	return(np.cos((np.pi*x)/2))


t,F,M = Mat(-1,1,K_test,F_test,10)
U1 = np.linalg.solve(M,F)
U =list()
for i in t:
	u=U_test(i)
	U.append(u)

error = np.linalg.norm(np.transpose(U1)-U)
print(error)
print("U=",U1)
plt.figure("Présentation de benchmark utilisé")
plt.plot(t,U1, label="solution approchée")
plt.plot(t,U, label="Solution exacte")
plt.legend()
plt.show()



#3.3

def K_electro_Love(x,t):
	return((1/np.pi)*(1/(1+(x-t)**2)))

def F_electro_Love(x):
	return(1)



t1,F1,M1 = Mat(-1,1,K_electro_Love,F_electro_Love,10)
U2 = np.linalg.solve(M1,F1)
print("U =",U2)
plt.figure("Equation Love en electrostatique")
plt.plot(t1,U2,label='Love en electrostatique')
plt.legend()
plt.show()


#4
h = 0.0002

def RLCprim(Y,t):
	e = 10
	R = 3
	L = 0.5
	C = 10E-6
	Y_p = np.zeros(2)
	Y_p[0] = (1/L)*(e- R*Y[0] - Y[1])
	Y_p[1] = (1/C)*Y[0]
	return Y_p

Y0 = np.array([0,0])
t2 = np.arange(0,2,0.0002)


Yodeint = sp.odeint(RLCprim, Y0,t2)

plt.figure("RLC circuit current intensity")
plt.plot(t2, Yodeint[:,0], color="hotpink", label='current intensity')
plt.xlabel('time (s)')
plt.legend()
plt.show()

plt.figure("RLC circuit voltage")
plt.plot(t2,Yodeint[:,1], color="sandybrown", label='voltage')
plt.xlabel('time (s)')
plt.legend()
plt.show()

#5 

def CCmotor(Y,t):
	R = 5
	L = 50E-3
	K_e = 0.2
	K_c = 0.1
	F_m =0.01 
	J_m = 0.05
	if (t>=10 and t<=50):
		u = 5
	else :
		u =0
	Y_p = np.zeros(2)
	Y_p[0] = (1/L)*(u - R*Y[0] - K_e*Y[1])
	Y_p[1] = (1/J_m)*(K_c*Y[0] - F_m*Y[1])
	return Y_p


Y_20 = np.array([0,0])
t3 = np.arange(0,80,0.002)

Yodeint_2 = sp.odeint(CCmotor,Y_20,t3)

plt.figure("RCC motor torque")
plt.plot(t3,Yodeint_2[:,0], color="mediumturquoise", label="torque")
plt.xlabel("time (s)")
plt.legend()
plt.show()

plt.figure("RCC motor angular velocity")
plt.plot(t3, Yodeint_2[:,1], color="darkorchid", label="angular velocity")
plt.xlabel("time (s)")
plt.legend()
plt.show()

#6

def rocket(Y,t):
	D = 4
	a = 8*10**3
	g = 9.81
	k =0.1
	u = 2*10**3

	Y_p = np.zeros(3)
	if (Y[1]<80):
		Y[1] = 80
		D =0.

	Y_p[0] =(D*u)/(Y[1]) - g - k*np.exp(-(Y[2])/(a))*(Y[0]**2)/(Y[1])
	Y_p[1] = -D 
	Y_p[2] = Y[0]
	return Y_p

Y_30 = np.array([0,400,0])
t4 = np.arange(0,160,0.1)
Yodeint_3 = sp.odeint(rocket,Y_30,t4)

plt.figure("Rocket pace")
plt.plot(t4, Yodeint_3[:,0], color="limegreen", label="Rocket pace")
plt.xlabel("time (s)")
plt.legend()
plt.show()

plt.figure("Rocket trajectory")
plt.plot(t4, Yodeint_3[:,2], color="darkorange", label="rocket trajectory")
plt.xlabel("time (s)")
plt.legend()
plt.show()

plt.figure("Rocket weight")
plt.plot(t4, Yodeint_3[:,1], color="deepskyblue", label="Rocket weight")
plt.xlabel("time (s)")
plt.legend()
plt.show()

#7 

def Predator_prey_model(Y,t):
	alpha_1 = 3
	beta_1 = 1
	alpha_2 = 2 
	beta_2 = 1

	Y_p = np.zeros(2)
	Y_p[0] = alpha_1*Y[0] - beta_1*Y[0]*Y[1]
	Y_p[1] = -alpha_2*Y[1] + beta_2*Y[0]*Y[1]

	return Y_p 


#Q1

t5 = np.arange(0,10,0.01)
Y1_0 = np.array([5,0])

Yp_odeint = sp.odeint(Predator_prey_model, Y1_0,t5)

plt.figure("Prey evolution w/out predator")
plt.plot(t5, Yp_odeint[:,0], color="springgreen", label="prey evolution w/out predator")
plt.plot(t5, Yp_odeint[:,1], color="dodgerblue", label="predator quantity")
plt.xlabel("time (s)")
plt.legend()
plt.show()


Y2_0 = np.array([0,3])

Yp2_odeint = sp.odeint(Predator_prey_model, Y2_0,t5)

plt.figure("Predator evolution w/out prey")
plt.plot(t5, Yp2_odeint[:,0], color="deeppink", label="prey quantity")
plt.plot(t5, Yp2_odeint[:,1], color="steelblue", label="predator evolution w/out prey")
plt.xlabel("time (s)")
plt.legend()
plt.show()


def Euler(Y0, t ,f):
	n= np.shape(Y0)
	N = len(t)
	Ye = np.zeros((N,n[0]))
	Ye [0,:] = Y0
	for i in range(N-1):
		Y = Ye[i,:]
		Y_prime = f(Y,t)
		Ye[i+1,:] = Y + h*Y_prime

	return Ye

Y3_0 = np.array([5,3])
h = 0.01
Ye = Euler(Y3_0, t5, Predator_prey_model)

plt.figure("Predator prey model w/ Euler")
plt.plot(t5, Ye[:,0], color="goldenrod", label="prey evolution w/ Euler explicit")
plt.plot(t5, Ye[:,1], color="mediumpurple", label="predator evolution w/ Euler explicit")
plt.xlabel("time (s)")
plt.legend()
plt.show()


Y4_0 = np.array([5,3])
Yp3_odeint = sp.odeint(Predator_prey_model, Y4_0, t5)

plt.figure("Predator prey model w/ odeint")
plt.plot(t5, Yp3_odeint[:,0], color="lightcoral", label="prey evolution w/ odeint")
plt.plot(t5, Yp3_odeint[:,1], color="slateblue", label="predator evolution w/ odeint")
plt.xlabel("time")
plt.legend()
plt.show()


plt.figure("Phase portrait")
plt.plot(Ye[:,0], Ye[:,1], color="tomato", label="Euler phase portrait")
plt.plot(Yp3_odeint[:,0], Yp3_odeint[:,1], color="blueviolet", label="Odeint phase portrait")
plt.xlabel("$y_1(t)$")
plt.legend()
plt.show()



#Code generalization for mutliples tests 

print("choisissez des valeurs pour alpha_1, beta_1 , alpha_2 , beta_2")
a = float(input("alpha_1 = "))
b = float(input("alpha_2 = "))
c = float(input("beta_1 = "))
d = float(input("beta_2 = "))

print("Choissisez des valeurs pour Y[0] et Y[1]")

p1 = float(input("Y[0] = "))
p2 = float(input("Y[1] = "))

def Predator_prey_model_2(Y,t):
	alpha_1 = a
	beta_1 = c
	alpha_2 = b
	beta_2 = d

	Y_p = np.zeros(2)
	Y_p[0] = alpha_1*Y[0] - beta_1*Y[0]*Y[1]
	Y_p[1] = -alpha_2*Y[1] + beta_2*Y[0]*Y[1]

	return Y_p 


Y5_0 = np.array([p1,p2])
Yp3_odeint = sp.odeint(Predator_prey_model_2, Y5_0, t5)

plt.figure("Predator prey model w/ odeint")
plt.plot(t5, Yp3_odeint[:,0], color="mediumblue", label="prey evolution w/ odeint")
plt.plot(t5, Yp3_odeint[:,1], color="magenta", label="predator evolution w/ odeint")
plt.xlabel("time")
plt.legend()
plt.show()


