#CS 6347 - Spring 2023
#Homework 2, Problem 2
#BY: JAYITA MALIK

import numpy as np

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# HELPER FUNCTIONS 

#Initialize the message 3D matrix
def initializeMessages(A, W):
    rows, cols = A.shape
    num_colors = W.shape[0]
    mess_arr = np.zeros((num_colors, rows, cols))

    for xj in range(num_colors):
        for i in range(rows):
            for j in range(cols):
                    if(A[i][j]==1):
                         mess_arr[xj][i][j] = 1
    return mess_arr

#Initialize beliefs matrix
def initializeBeliefs():
    b_arr = np.zeros((A1.shape[0], W1.shape[0]))
    return b_arr

#Get neighbors of a vertex
def neigh(i):
    my_list = []
    for m in range(A1.shape[0]):
        if(A1[i][m]==1):
            my_list.append(m)
    return my_list

#Get neighbors of a vertex excluding j
def neigh_except_j(i, j):
    my_list = []
    for m in range(A1.shape[0]):
        if(A1[i][m]==1 and m!=j):
            my_list.append(m)
    return my_list

#Get the product of incoming messages to a vertex
def getProd(i, xi, messages):
    prod = 1
    neighbors_of_i = neigh(i)
    for k in neighbors_of_i:
        prod = prod * messages[xi][k][i]
    return prod

#Get the product of incoming messages (except from j)
def getProd_except_j(i, j, messages, xi):
    prod = 1
    neighbors_of_i_except_j = neigh_except_j(i, j)
    for k in neighbors_of_i_except_j:
        prod = prod * messages[xi][k][i]
    return prod

#Get psi
def getPsi(xi, xj):
    if(xi==xj):
        return 0
    else:
        return 1
    
#Get phi
def getPhi(xi):
    x = np.exp(W1[xi])
    return x

#Get scaling factor for messages
def getScalingFactor(W, messages, i, j):
    running_sum = 0
    for xj in range(W.shape[0]):
        for xi in range(W.shape[0]):
            running_sum += (getProd_except_j(i, j, messages, xi) * getPsi(xi, xj) * getPhi(xi))
    sf = running_sum
    return sf

#Create a new message for sum-product
def getNewMessage_SP(W, i, j, messages, xj):
    sum = 0
    for xi in range(W.shape[0]):
        prod = getProd_except_j(i, j, messages, xi)
        psi = getPsi(xi, xj) 
        phi = getPhi(xi)
        total_prod = prod * psi * phi
        sum = sum + total_prod
    return sum

#Create a new message for max-product
def getNewMessage_MP(W, i, j, messages, xj):
    a = []
    for xi in range(W.shape[0]):
        phi_i = getPhi(xi)
        psi_ij = getPsi(xi, xj)
        prod = getProd_except_j(i, j, messages, xi)
        total_prod = phi_i * psi_ij * prod
        a.append(total_prod)
    max_xi = np.max(np.array(a))
    return max_xi

#Update the 3D message matrix
def updateMessage(i, j, xj, messages, new_message):
    messages[xj][i][j] = new_message

#Get normalization constant 
def getNormalConst(messages, i):
    running_sum = 0
    for xi in range(W1.shape[0]):
            running_sum += (getPhi(xi)*getProd(i, xi, messages))
    nc = running_sum
    return nc

#Returns belief of a vertex i
def getBelief(i, messages):
    b_arr = np.zeros((W1.shape[0]))
    normal_const = getNormalConst(messages, i)
    for xi in range(W1.shape[0]):
        phi = getPhi(xi)
        prod = getProd(i, xi, messages)
        bi_xi = (phi * prod)/normal_const
        b_arr[xi] = bi_xi
    return b_arr

#Returns the coloring assignment
def getColoringAssignment(beliefs):
    asgnmt_list = []
    for belief in beliefs:
        x = np.argwhere(belief == np.amax(belief))
        x = x.flatten().tolist()
        if(len(x)>1):
            color = 0
        else:
            color = x[0]+1
        asgnmt_list.append(color)
    return asgnmt_list

#Gets the entropy for bethe free energy
def getEntropy(i, beliefs):
    summation = 0
    for xi in range(W1.shape[0]):
        belief_xi = beliefs[i][xi]
        if(belief_xi==0):
            log_belief_xi = 0
        else:
            log_belief_xi = np.log(belief_xi)
        summation += (belief_xi * log_belief_xi)
    return -summation

#Gets normalization constant for pairwise beliefs
def getNormalConst_pw(i, j, messages):
    running_sum = 0
    for xi in range(W1.shape[0]):
        for xj in range(W1.shape[0]):
            running_sum += (getPhi(xi)*getPhi(xj)*getPsi(xi,xj)*
                            getProd_except_j(i, j, messages, xi)*getProd_except_j(j, i, messages, xj))
    nc = running_sum
    return nc

#Calculates pairwise beliefs dynamically 
def getPairwiseBelief(i, j, xi, xj, messages):
    phi_i = getPhi(xi)
    phi_j = getPhi(xj)
    psi_ij = getPsi(xi, xj)
    prod_without_j = getProd_except_j(i, j, messages, xi)
    prod_without_i = getProd_except_j(j, i, messages, xj)
    pw_belief = phi_i * phi_j * psi_ij * prod_without_i * prod_without_j

    nc = getNormalConst_pw(i, j, messages)
    normalized_pw_belief = pw_belief/nc

    return normalized_pw_belief

#Returns the second term (I_ij) for Bethe Free Energy
def getEnergy(i, j, messages, beliefs):
    energy_sum = 0
    for xi in range(W1.shape[0]):
        for xj in range(W1.shape[0]):
            if(xi==xj):
                continue
            if(i==j or A1[i][j]==0):
                pw_belief_xi_xj = 0
            else:
                pw_belief_xi_xj = getPairwiseBelief(i, j, xi, xj, messages)
            if(pw_belief_xi_xj==0):
                log_pw_belief = 0
            else:
                fraction = pw_belief_xi_xj/(beliefs[i][xi]*beliefs[j][xj])
                log_pw_belief = np.log(fraction)
            energy_sum = energy_sum + (pw_belief_xi_xj * log_pw_belief)
    return energy_sum

#Gets Bethe Free Energy
def bethe_free_energy(beliefs, messages):
    big_entropy = 0
    for i in range(A1.shape[0]):
        entropy_i = getEntropy(i, beliefs)
        big_entropy += entropy_i
    
    big_energy_sum = 0 
    for i in range(A1.shape[0]):
        for j in range(i+1, A1.shape[0]):
            energy_ij = getEnergy(i, j, messages, beliefs)
            big_energy_sum += energy_ij

    bfe = big_entropy-big_energy_sum
    return bfe

#Gets partition function Z
def getPartitionFunc(beliefs, messages):
    log_z = bethe_free_energy(beliefs, messages)
    z = np.exp(log_z)
    return z

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# SUM PRODUCT FUNCTION
def sumprod(A, W, its):
    r, c = A.shape
    messages_t = initializeMessages(A, W)
    messages_t_minus_1 = initializeMessages(A, W)
    beliefs_t = initializeBeliefs()

    for t in range(its):
        for i in range(r):
            for j in range(c):
                if(A1[i][j]==1):
                    scaling_factor = getScalingFactor(W, messages_t, i, j)
                    for xj in range(W.shape[0]):
                        new_message_val = getNewMessage_SP(W, i, j, messages_t_minus_1, xj)
                        scaled_message = new_message_val/scaling_factor
                        updateMessage(i, j, xj, messages_t, scaled_message) 
        for i in range(r):
            beliefs__for_i = getBelief(i, messages_t)
            beliefs_t[i] = beliefs__for_i

        messages_t_minus_1 = np.copy(messages_t)
    
    assignment_list = getColoringAssignment(beliefs_t)
    print("Coloring assignment: ", assignment_list)

    Z = getPartitionFunc(beliefs_t, messages_t)
    return Z
#----------------------------------------------------------------------------------


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MAX PRODUCT FUNCTION                   
def maxprod(A, W, its):
    r, c = A.shape
    messages_t = initializeMessages(A, W)
    messages_t_minus_1 = initializeMessages(A, W)
    beliefs_t = initializeBeliefs()
    for t in range(its):
        for i in range(r):
            for j in range(c):
                if(A[i][j]==1):
                    for xj in range(W.shape[0]):
                        new_message_val = getNewMessage_MP(W, i, j, messages_t_minus_1, xj)
                        updateMessage(i, j, xj, messages_t, new_message_val)
        for i in range(r):
            beliefs__for_i = getBelief(i, messages_t)
            beliefs_t[i] = beliefs__for_i

        messages_t_minus_1 = np.copy(messages_t)
    
    assignment_list = getColoringAssignment(beliefs_t)
    return assignment_list
#----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# MAIN
global A1
global W1
print("\n------------------------------")
print("Would you like to use an existing example, or enter your own matrix, weights, and iterations?")
source = input("e - Example\no - Own matrix, weights, iterations\n")

#If user wants to use the in-built examples
if(source=="E" or source=="e"):

    #Example 1
    print("\n~~***----  EXAMPLE 1  ----***~~")
    A1 = np.array([[0,1,0], 
                   [1,0,1], 
                   [0,1,0]])
    W1 = np.array([1,2])
    its1 = 30

    print("GRAPH:   1--2--3")
    print("WEIGHTS: [1,2]")
    print("ITS:     30\n")

    print("\n------- SUM PRODUCT --------")
    Z = sumprod(A1, W1, its1)
    print("Approximate parition function, Z = ", Z)
    print("----------------------------\n")

    print("------- MAX PRODUCT --------")
    assignment_list = maxprod(A1, W1, its1)
    print("Coloring assignment for [v1,v2,...,vn]: ", assignment_list)
    print("----------------------------\n")

    #Ask user if they want to see 1 more example
    more_example = input("Would you like to see 1 more example? (y/n)\n")
    if(more_example=="y" or more_example=="Y"):

        #Example 2
        print("\n~~***----  EXAMPLE 2  ----***~~")
        A1 = np.array([[0,1,0,0], 
                       [1,0,1,1], 
                       [0,1,0,0],
                       [0,1,0,0]])
        W1 = np.array([1,1])
        its1 = 30

        print("GRAPH:   1--2--3")
        print("            |   ")
        print("            4   ")
        print("WEIGHTS: [1,1]")
        print("ITS:     30\n")

        print("\n------- SUM PRODUCT --------")
        Z = sumprod(A1, W1, its1)
        print("Approximate parition function, Z = ", Z)
        print("----------------------------\n")

        print("------- MAX PRODUCT --------")
        assignment_list = maxprod(A1, W1, its1)
        print("Coloring assignment for [v1,v2,...,vn]: ", assignment_list)
        print("----------------------------\n")

#If the user wants to enter own matrix, weights, iter
elif(source=="O" or source=="o"):
    print()

    #Input matrix
    R = int(input("For matrix, enter the # of rows/cols: "))
    C = R
    M1 = []
    print("Enter the entries rowwise (press ENTER after each entry):")
    
    for i in range(R):          
        a =[]
        for j in range(C):      
            a.append(int(input()))
        M1.append(a)
    
    print("Your matrix is: ")
    for i in range(R):
        for j in range(C):
            print(M1[i][j], end = " ")
        print()
    A1 = np.array(M1)


    #Input weights
    v = int(input("\nFor weights, enter the # of values:"))
    M2 = []

    print("Enter the weights (press ENTER after each entry):")
    for i in range(v):
        M2.append(int(input()))

    print("Your weights are: ", M2)
    W1 = np.array(M2)


    #Input iterations
    its2=int(input("\nEnter # of iterations: "))


    #Sum-product and max-product functions
    print("\n------- SUM PRODUCT --------")
    Z = sumprod(A1, W1, its2)
    print("Approximate parition function, Z = ", Z)
    print("----------------------------\n")

    print("------- MAX PRODUCT --------")
    assignment_list = maxprod(A1, W1, its2)
    print("Coloring assignment for [v1,v2,...,vn]: ", assignment_list)
    print("----------------------------\n")


else:
    print("Exiting program: Invalid input, enter E or O")

