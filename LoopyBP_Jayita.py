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

def neigh_except_j(i, j):
    my_list = []
    for m in range(A1.shape[0]):
        if(A1[i][m]==1 and m!=j):
            my_list.append(m)
    #print("Neighbours of i=", i, " for given j=", j, " are: ", my_list)
    return my_list


def getProd(i, j, messages, xi):
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
        prod = getProd(i, j, messages, xi)
        psi = getPsi(xi, xj) 
        phi = getPhi(xi)
        total_prod = prod * psi * phi
        sum = sum + total_prod
    return sum

def updateMessage(i, j, xj, messages, new_message):
    print("!!Message update!!")
    messages[xj][i][j] = new_message
    messages[xj][j][i] = new_message #maybe problem??????????
    print("New matrix: ")
    print(messages)
    
def getScalingFactor(xj, messages):
    running_sum = 0
    rows, cols = A1.shape
    for i in range(rows):
        for j in range(cols):
            running_sum += messages[xj][i][j]
    sf=1/running_sum
    return sf

def sum_product(A, W, its):
    r, c = A.shape
    messages_t = initializeMessages(A, W)
    print("Initial message matrix: \n", messages_t)
    #Do we need a copy of messages_t?????????????????????????????

    for t in range(its):
         print("ITERATION -- ", t, " ----------------")
         for i in range(r):
              for j in range(c):
                for xj in range(W.shape[0]):
                    new_message_val = getNewMessage_SP(W, i, j, messages_t, xj)
                    #scaling_factor = 0.4 #PROBLEM???????????????????????????????????????
                    scaling_factor = getScalingFactor(xj, messages_t)
                    scaled_message = scaling_factor * new_message_val
                    updateMessage(i, j, xj, messages_t, scaled_message)
              


#-----------------------------------------------------------------------------------
#Graph 1: 1--2--3
#Weights: [1, 2]

global A1 
global W1 
A1 = np.array([[0,1,0], [1,0,1], [0,1,0]])
W1 = np.array([1,2])
its1 = 5

print("Your matrix: \n", A1)
sum_product(A1, W1, its1)
