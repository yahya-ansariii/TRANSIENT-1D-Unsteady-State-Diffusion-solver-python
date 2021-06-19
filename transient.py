'''      Copyright 2021 MOHAMMED YAHYA ANSARI

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
'''
import pandas as pd
import numpy as np


print("\n \t\t TRANSIENT Solver\n")
n = int(input("\n\tEnter the no. of grid points:   "))

l = float(input("\n\tEnter length of plate in m:   "))

tk = float(input("\n\tEnter thermal conductivity in W/mK or W/mC:   "))

rc = float(input("\n\tEnter rhoC in J/m3K:   "))

ti = float(input("\n\tEnter Initial temperature at T = 0 sec:   "))

te = float(input("\n\tEnter temperature at east side of plate:   "))

dt = float(input("\n\tEnter dt:   "))

step = int(input("\n\tEnter number of time step:   "))


# create empty list
D = [0]*n
beta = [0]*n
alpha = [0]*n
c = [0]*n
A = [0]*n
C = [0]*n
X = [0]*n
temp = [0]*n
tempold = [0]*n

dx = l/n
Z = (rc*dx)/dt
z = tk/dx

#TDMA function very specific to this numerical solver, for general tdma solver refer : https://github.com/novus-afk/TDMA-Solver
def TDMA(n, beta, D, alpha, c):
    beta[0] = 0
    beta[n-1] = beta[1]
    alpha[0] = alpha[1]
    alpha[n-1] = 0
    #copy common values
    for i in range(1, n-1):
        D[i] = D[1]
        beta[i] = beta[1]
        alpha[i] = alpha[1]
    # solve forward substitution
    for i in range(0, n):
        A[i] = alpha[i]/(D[i] - beta[i]*A[i-1])
        C[i] = (beta[i]*C[i-1] + c[i])/(D[i] - beta[i]*A[i-1])

    X[n-1] = C[n-1]
    # solve backward substitution
    j = n-2
    while j >= 0:
        X[j] = A[j] * X[j+1] + C[j]
        j = j-1

    return X


# create dataframe to store result in table
df = pd.DataFrame(columns=["1", "2", "3", "4", "5"])

# Switch case for type of numerical
choice = ""
while choice != "q":
    print("""\n\t\tSelect a Scheme :
        
        [ 1 ] Explicit Scheme

        [ 2 ] Implicit Scheme

        [ 3 ] Crank-Nicolson Scheme 

        [ q ] Exit\n""")
    choice = input("\n\tEnter Choice :\t")

    if choice == "1":
        print("\n\tExplicit Scheme\n")
        for o in range(0, n):
            temp[o] = ti
            tempold[o] = ti

        for i in range(0, step+1):
            print("Temperature values for time step {} are:\n" .format(i))
            for k in range(0, n):
                print(" T%d = %f C\t" % (k+1, temp[k]))
            print("\n")
            df.loc[len(df.index)] = tempold
            temp[0] = (((Z-z)*tempold[0])+(z*tempold[1]))/Z

            for j in range(1, n-1):
                temp[j] = (((Z-(2*z))*tempold[j]) +
                           (z*tempold[j-1])+(z*tempold[j+1]))/Z

            temp[n-1] = (((Z-(3*z))*tempold[n-1]) +
                         (z*tempold[n-2]) + (2*te*z))/Z
            tempold = temp.copy()
        break


    elif choice == "2":
        print("\n\tImplicit Scheme\n")
        #initialize data for tdma function
        D[0] = Z+z
        D[1] = Z+z+z
        D[n-1] = Z+z+z+z
        beta[1] = z
        alpha[1] = z
        #copy initial temperature
        for o in range(0, n):
            temp[o] = ti
        tempold = temp.copy()

        print("Temperature values are:\n")

        #print temp values
        for p in range(0, step+1):
            for i in range(0, n):
                print("T%d = %f C \t" % (i+1, temp[i]))
            print("\n")
            
            #add new temp values to the pandas dataframe(add new temperature row, similar to printing values at each time step)
            df.loc[len(df.index)] = temp

            for i in range(0, n):
                c[i] = Z*tempold[i]

            temp = TDMA(n, beta, D, alpha, c)
            tempold = temp.copy()
        break

    #solve for crank nicolson scheme (currently it is a copy of implicit)
    elif choice == "3":
        print("\n\tCrank Nicolson Scheme\n")
        D[0] = Z+(z/2)
        D[1] = Z+z
        D[n-1] = Z+z
        beta[1] = z/2
        alpha[1] = z/2
        
        for o in range(0, n):
            temp[o] = ti
        tempold = temp.copy()

        print("Temperature values are:\n")

        for p in range(0, step+1):
            for i in range(0, n):
                print("T%d = %f C \t" % (i+1, temp[i]))
            print("\n")
            df.loc[len(df.index)] = temp

            #set constant values
            c[0] = (Z * tempold[0]) + (z * 0.5 * (tempold[1] - tempold[0]))
            c[n-1] = ((Z-z)*tempold[n-1]) + (z * 0.5 * tempold[n-2])
            for i in range(1, n-1):
                c[i] = ((Z-z)*tempold[i]) + (z * 0.5 * (tempold[i-1]+tempold[i+1]))

            temp = TDMA(n, beta, D, alpha, c)
            tempold = temp.copy()
        break

    elif choice == "q":
        exit()

    else:
        print("\n\n\tInvalid choice, Try again!\n")


print(df)

input("Press enter to exit")
