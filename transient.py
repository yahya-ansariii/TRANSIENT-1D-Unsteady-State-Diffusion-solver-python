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

# take user input
print("\n \t\t TRANSIENT Solver\n")
print("""\n\tTo learn more about Transient heat conduction read the Transient.pdf
    
    (The left boundary is insulated in this solver)
    """)
n = int(input("\n\tEnter the no. of grid points :   "))

l = float(input("\n\tEnter length of plate in m :   "))

tk = float(input("\n\tEnter thermal conductivity in W/mK or W/mC :   "))

rc = float(input("\n\tEnter \N{GREEK SMALL LETTER RHO}C in J/m3K :   "))

ti = float(input("\n\tEnter initial temperature at T = 0 sec :   "))

te = float(input("\n\tEnter temperature at east side of plate :   "))

dt = float(input("\n\tEnter dt :   "))

step = int(input("\n\tEnter number of time step :   "))


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

# create dataframe to store result in table
df = pd.DataFrame(columns=["Node: 1", "Node: 2",
                  "Node: 3", "Node: 4", "Node: 5"])

dx = l/n
Z = (rc*dx)/dt
z = tk/dx

# Define TDMA function very specific to this numerical solver, for general tdma solver refer : https://github.com/novus-afk/TDMA-Solver


def TDMA(n, beta, D, alpha, c):
    beta[0] = 0
    beta[n-1] = beta[1]
    alpha[0] = alpha[1]
    alpha[n-1] = 0
    # copy common values
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


def EXPLICIT(temp, tempold):
    for o in range(0, n):
        temp[o] = ti
        tempold[o] = ti

    for i in range(0, step+1):
        print("Temperature values for time step {} are:\n" .format(i))
        for k in range(0, n):
            print(" T%d = %f C\t" % (k+1, temp[k]))
        print("\n")
        df.loc[len(df.index)] = tempold  # add data to dataframe last row
        temp[0] = (((Z-z)*tempold[0])+(z*tempold[1]))/Z

        for j in range(1, n-1):
            temp[j] = (((Z-(2*z))*tempold[j]) +
                       (z*tempold[j-1])+(z*tempold[j+1]))/Z

        temp[n-1] = (((Z-(3*z))*tempold[n-1]) + (z*tempold[n-2]) + (2*te*z))/Z
        tempold = temp.copy()


# noinspection PyUnresolvedReferences
def IMPLICIT(temp, tempold):
    # initialize data for tdma function
    D[0] = Z+z
    D[1] = Z+z+z
    D[n-1] = Z+z+z+z
    beta[1] = z
    alpha[1] = z
    # copy initial temperature
    for o in range(0, n):
        temp[o] = ti
    tempold = temp.copy()

    for p in range(0, step+1):
        print("Temperature values for time step {} are:\n" .format(p))
        # print temp values
        for i in range(0, n):
            print("T%d = %f C \t" % (i+1, temp[i]))
        print("\n")

        # add new temp values to the pandas dataframe(add new temperature row, similar to printing values at each time step)
        df.loc[len(df.index)] = temp

        for i in range(0, n):
            c[i] = Z*tempold[i]
        temp = TDMA(n, beta, D, alpha, c)
        tempold = temp.copy()


def CrankN(temp, tempold):
    D[0] = Z+(z/2)
    D[1] = Z+z
    D[n-1] = Z+z+z/2
    beta[1] = z/2
    alpha[1] = z/2

    for o in range(0, n):
        temp[o] = ti
    tempold = temp.copy()

    print("Temperature values are:\n")

    for p in range(0, step+1):
        print("Temperature values for time step {} are:\n" .format(p))
        # print temp values
        for i in range(0, n):
            print("T%d = %f C \t" % (i+1, temp[i]))
        print("\n")
        # add result to DataFrame
        df.loc[len(df.index)] = temp

        # set constant values
        c[0] = (Z * tempold[0]) + (z * 0.5 * (tempold[1] - tempold[0]))
        c[n-1] = ((Z-z-z/2)*tempold[n-1]) + (z * tempold[n-2]) + (2*z * te)
        for i in range(1, n-1):
            c[i] = ((Z-z)*tempold[i]) + (z * 0.5 * (tempold[i-1]+tempold[i+1]))

        temp = TDMA(n, beta, D, alpha, c)
        tempold = temp.copy()


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
        EXPLICIT(temp, tempold)
        break

    elif choice == "2":
        print("\n\tImplicit Scheme\n")
        IMPLICIT(temp, tempold)
        break

    # solve for crank nicolson scheme (currently it is a copy of implicit)
    elif choice == "3":
        print("\n\tCrank Nicolson Scheme\n")
        CrankN(temp, tempold)
        break

    elif choice == "q":
        exit()

    else:
        print("\n\n\tInvalid choice, Try again!\n")


# add time column
Ts = []
for i in range(0, step+1):
    Ts.append(i*dt)
df.insert(0, "Time(sec)", Ts)

# print DataFrame
print(df)


# save choice
while choice != "q":
    print('''\n\t[ y ] Save result to excel file.
    
        [ q ] Exit without saving result
        ''')
    choice = input("\nEnter yout choice :\t")
    if choice == "y":
        # add time step column at the start of the DataFrame
        df.insert(0, 'Time Step.', range(0, len(df)))
        df.to_excel("output//Transient.xlsx", sheet_name='Output',
                    index=False)  # .to_excel to export excel file
        print(
            "\n\n*************** Export complete! Check output folder. ***************\n\n")
        break
    elif choice == "q":
        print("\n***** Result not saved *****")
        break
    else:
        print("\n Invalid Choice, Try again !")


input("\nPress Enter to exit")
