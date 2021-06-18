import pandas as pd
import numpy as np


print("\n \t\t TRANSIENT Solver\n")
n = int(input("\nEnter the no. of grid points: "))

l = float(input("Enter length of plate in m: "))

tk = float(input("Enter thermal conductivity in W/mK or W/mC: "))

rc = float(input("Enter rhoC in J/m3K: "))

ti = float(input("Enter Initial temperature at T = 0 sec: "))

te = float(input("Enter temperature at east side of plate: "))

dt = float(input("Enter dt: "))

step = int(input("Enter number of time step: "))





#create empty list
D = [0]*n
beta = [0]*n
alpha = [0]*n
c =[0]*n
A = [0]*n
C = [0]*n
temp = [0]*n
tempold = [0]*n
dx = l/n
Z= (rc*dx)/dt
z= tk/dx

print("\n\t\t\t**MENU**")
print("\n\t1) EXPLICIT")
print("\n\t2) IMPLICIT")
print("\n\t3) EXIT")
cho = int(input("\n\t ENTER YOUR CHOICE: "))


df = pd.DataFrame(columns = ["1", "2", "3", "4" ,"5"])

#switch using if and elseif statements

if cho == 1:
    for o in range(0 , n):
        temp[o]=ti
        tempold[o]=ti
    
    for i in range(0 , step+1):
        print("Temperature values for time step {} are:\n" .format(i))
        for k in range(0 , n):
            print(" T%d = %f C\t" % (k+1, temp[k]))
        print("\n")
        df.loc[len(df.index)] = tempold 
        temp[0] = (((Z-z)*tempold[0])+(z*tempold[1]))/Z

        for j in range(1 , n-1):
            temp[j]= (((Z-(2*z))*tempold[j])+(z*tempold[j-1])+(z*tempold[j+1]))/Z

        temp[n-1]= (((Z-(3*z))*tempold[n-1])+(z*tempold[n-2])+ (2*te*z))/Z
        tempold = temp.copy()


    



elif cho == 2:
    D[0]=Z+z
    D[1]=Z+z+z
    D[n-1]=Z+z+z+z
    beta[1]=z
    alpha[1]=z
    for o in range(0 , n):
        temp[o]=ti

    print("Temperature values are:\n")

    for p in range(0 , step+1):
        for i in range(0 , n):
            print("T%d = %f C \t" % (i+1, temp[i]))
        print("\n")
        df.loc[len(df.index)] = temp

        for i in range(0 , n):
            c[i] = Z*temp[i]

        beta[0] = 0
        beta[n-1] = beta[1]
        alpha[0] = alpha[1]
        alpha[n-1] = 0

        for i in range(2 , n-1):
            D[i] = D[1]
            beta[i]=beta[1]
            alpha[i]=alpha[1]



        a = alpha[0]
        A[0] = a/D[0]
        C[0] = c[0]/D[0]
        for i in range(1 , n):
            A[i] = alpha[i]/(D[i] - beta[i]*A[i-1])
            C[i] = (beta[i]*C[i-1] + c[i])/(D[i] - beta[i]*A[i-1])
    
        temp[n-1] = C[n-1]

        j = n-2
        while j>=0:
            temp[j] = A[j] * temp[j+1] + C[j]
            j = j-1

elif cho == 3:
    print("Closing program")

else:
    print("invalid option")



print(df)
input("Press enter to exit")