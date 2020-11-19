
# Alex Tatievsky 317414506 #
# https://github.com/Alextatie/numericalHW

from numpy import array, zeros, diag, diagflat, dot, sum, all, shape, linalg,ones


def diagcheck(A):
    temp=A
    D = diag(abs(temp))  # diagonal numbers in abs value
    S = sum(abs(temp), axis=1) - D  # sum of rows without the diagonal number in abs value
    if all(D >= S):  # checks in every row if the diagonal number is >= the rest in abs value
        return(True)
    else:
        return(False)

def swapcheck(A,b):
    dim = shape(A)[0]
    keys = list(-1 for _ in range(dim))  # making an empty list to containt row numbers
    for i in range(dim):
        for j in range(dim):
            if a[i][j] >= sum(a[i]) - a[i][j]:
                if keys[i] == -1:
                    keys[i] = j
    if len(keys) > len(set(keys)):  # checking if every item is unique, if so, able to swap rows
        return (False,False)
    else:
        newA = zeros(shape(A))
        newb = zeros(shape(b))
        i = 0
        for key in keys:
            newA[key] = A[i]
            newb[key] = b[i]
            i = i + 1
        return(newA, newb)

def LUD(A):
    dim=shape(A)[0]  #dimension via built in func
    D = diagflat(diag(a))  #diagonal via built in funcs
    L = zeros(shape(A), dtype=float)  #a loop to create L by copying the value from A if i>j (meaning all the lower slots)
    for i in range (dim):
        for j in range (dim):
            if i>j:
                L[i][j]=A[i][j]
    U = A-L-D
    return(L,U,D)

def Jacobi(A,b):
    L = A[0]
    U = A[1]
    D = A[2]
    dim = shape(L)[0]
    x = zeros([dim, 1]) #initiation with 0's
    xa = ones([dim, 1]) #initiation of Xr-1 for the condition
    i = 0 # index
    s = 0 # stop flag
    print("X(", end='')
    print(i,end=') = (')
    for d in range(dim):
        if d==dim-1:
            print(round(float(x[d]),4),end=")\n")
        else:
            print(round(float(x[d]),4),end=', ')
    while s!=1:
        i = i + 1
        xa = x
        x = -dot(linalg.inv(D), dot((L + U), x)) + dot(linalg.inv(D), b)  # Jacobi formula
        s = 1
        print("X(", end='')    #long print function just so it will look better
        print(i, end=') = (')
        for d in range(dim):
            if abs(x[d] - xa[d]) > 0.001:
                s = 0
            if d == dim - 1:
                print(round(float(x[d]),4), end=")\n")
            else:
                print(round(float(x[d]),4), end=', ')
        if s == 1:
            break  #breaks if we dont find a difference larger than 0.001

def GaussSeidel(A,b):
    L = A[0]
    U = A[1]
    D = A[2]
    dim = shape(L)[0]
    x = zeros([dim, 1]) #initiation with 0's
    xa = ones([dim, 1]) #initiation of Xr-1 for the condition
    i = 0 # index
    s = 0 # stop flag
    print("X(", end='')
    print(i,end=') = (')
    for d in range(dim):
        if d==dim-1:
            print(round(float(x[d]),4),end=")\n")
        else:
            print(round(float(x[d]),4),end=', ')
    while s!=1:
        i = i + 1
        xa = x
        x= -dot(linalg.inv(L+D),dot((U),x))+dot(linalg.inv(L+D),b)  #Gauss-Seidel formula
        s = 1
        print("X(", end='')    #long print function just so it will look better
        print(i, end=') = (')
        for d in range(dim):
            if abs(x[d] - xa[d]) > 0.001:
                s = 0
            if d == dim - 1:
                print(round(float(x[d]),4), end=")\n")
            else:
                print(round(float(x[d]),4), end=', ')
        if s == 1:
            break  #breaks if we dont find a difference larger than 0.001

# main
dim=int(input("Enter a dimension for a quadratic matrix:\n"))
a=zeros([dim,dim])
print("Enter values row by row:\n")
for i in range(dim):
    for j in range(dim):
        a[i][j]=float(input())
b = zeros([dim, 1])
print("Enter values for matrix b:\n")
for i in range(dim):
    b[i][0] = float(input())
if diagcheck(a) or (not diagcheck(a) and type(swapcheck(a,b)[0]) == type(a)):
    print("Your matrix is:")
    for i in range(dim):
        print("(", end='')
        for j in range(dim):
            print(a[i][j], end=' ')
        print("|", b[i][0], ")")
    if not diagcheck(a):   #if we swapped, change our matrix's into the new ones that we swapped into
        s = swapcheck(a,b)
        a = s[0]
        b = s[1]
        print("\nYour rearanged matrix is:")
        for i in range(dim):
            print("(", end='')
            for j in range(dim):
                print(a[i][j], end=' ')
            print("|", b[i][0], ")")
    o = input("\n\nChoose a method:\n(a) Jacobi\n(b) GaussSeidel\n") # method choice
    if o == 'a':
        Jacobi(LUD(a), b)
    elif o == 'b':
        GaussSeidel(LUD(a), b)
    else:
        print("Wrong input.\n")
else:
    print("Swaping rows cant make a diagonally dominant matrix")
