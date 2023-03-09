import numpy as np

def initializeMessages(A, W):
    rows, cols = A.shape
    num_colors = W.shape[0]
    mess_arr = np.zeros((num_colors, rows, cols))

    #Make every edge as 1 
    for xj in range(num_colors):
        for i in range(rows):
            for j in range(cols):
                    if(A[i][j]==1):
                         mess_arr[xj][i][j] = 1
                         #print("Made message = 1 for i = ", i, " j = ", j)
                         #print("i = ", i, " j = ", j, " xj = ", xj, " || MES = ", mess_arr[xj][i][j])
    return mess_arr

def initializeBeliefs():
    b_arr = np.zeros((A1.shape[0], W1.shape[0]))
    return b_arr

def neigh(i):
    my_list = []
    for m in range(A1.shape[0]):
        if(A1[i][m]==1):
            my_list.append(m)
    return my_list

def neigh_except_j(i, j):
    my_list = []
    for m in range(A1.shape[0]):
        if(A1[i][m]==1 and m!=j):
            my_list.append(m)
    #print("Neighbours of i=", i, " for given j=", j, " are: ", my_list)
    return my_list

def getProd(i, xi, messages):
    prod = 1
    neighbors_of_i = neigh(i)
    for k in neighbors_of_i:
        prod = prod * messages[xi][k][i]
    return prod

def getProd_except_j(i, j, messages, xi):
    prod = 1
    neighbors_of_i_except_j = neigh_except_j(i, j)
    for k in neighbors_of_i_except_j:
        prod = prod * messages[xi][k][i]
    return prod

def getPsi(xi, xj):
    if(xi==xj):
        return 0
    else:
        return 1
    
def getPhi(xi):
    x = np.exp(W1[xi])
    return x

def getNewMessage_SP(W, i, j, messages, xj):
    sum = 0
    for xi in range(W.shape[0]):
        prod = getProd_except_j(i, j, messages, xi)
        psi = getPsi(xi, xj) 
        phi = getPhi(xi)
        total_prod = prod * psi * phi
        sum = sum + total_prod
    return sum

def updateMessage(i, j, xj, messages, new_message):
    #print("!!Message update!!")
    messages[xj][i][j] = new_message
    #messages[xj][j][i] = new_message #maybe problem??????????
    #print("New matrix: ")
    #print(messages)
    
def getScalingFactor(xj, messages):
    running_sum = 0
    rows, cols = A1.shape
    for i in range(rows):
        for j in range(cols):
            running_sum += messages[xj][i][j]
    sf=1/running_sum
    return sf

def getNormalConst():
    return 0.4

def getBelief(i, messages):
    b_arr = np.zeros((W1.shape[0]))
    normal_const = getNormalConst()
    for xi in range(W1.shape[0]):
        phi = getPhi(xi)
        print("Phi   == ", phi)
        prod = getProd(i, xi, messages)
        print("Prod  == ", prod)
        bi_xi = normal_const * phi * prod
        print("Bi_Xi == ", bi_xi)
        print()
        b_arr[xi] = bi_xi
    return b_arr



def sum_product(A, W, its):
    r, c = A.shape
    messages_t = initializeMessages(A, W)
    beliefs_t = initializeBeliefs()
    #print("Initial message matrix: \n", messages_t)
    #Do we need a copy of messages_t????????????????????????????? messages_t+1

    for t in range(its):
        print("ITERATION -- ", t, " ----------------")
        for i in range(r):
            for j in range(c):
                if(A1[i][j]==1):
                    for xj in range(W.shape[0]):
                        new_message_val = getNewMessage_SP(W, i, j, messages_t, xj)
                        scaling_factor = 0.4 #PROBLEM???????????????????????????????????????
                        #scaling_factor = getScalingFactor(xj, messages_t)
                        #print("SCALING FACTOR === ", scaling_factor)
                        scaled_message = scaling_factor * new_message_val
                        updateMessage(i, j, xj, messages_t, scaled_message) #REVIEW THIS
                        #xi = xj
                        #belief_xi = getBelief(xi)
            beliefs__for_i = getBelief(i, messages_t)
            beliefs_t[i] = beliefs__for_i
    
    print("Belief matrix becomes: \n")
    print(beliefs_t)
                        




#-----------------------------------------------------------------------------------
#Graph 1: 1--2--3
#Weights: [1, 2]

global A1 
global W1 
# A1 = np.array([[0,1,0], [1,0,1], [0,1,0]])
# W1 = np.array([1,2])
A1 = np.array([
    [0, 1, 0, 0],
    [1, 0, 1, 0],
    [0, 1, 0, 1],
    [0, 0, 1, 0]
])
W1 = np.array([1, 1, 1])
its1 = 10

#print("Your matrix: \n", A1)
sum_product(A1, W1, its1)
