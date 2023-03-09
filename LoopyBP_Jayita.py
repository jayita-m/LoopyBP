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

# def initializePwBeliefs():
    

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

def getScalingFactor(W, messages, i, j):
    running_sum = 0
    for xj in range(W.shape[0]):
        for xi in range(W.shape[0]):
            running_sum += (getProd_except_j(i, j, messages, xi) * getPsi(xi, xj) * getPhi(xi))
    sf = running_sum
    return sf

def getNewMessage_SP(W, i, j, messages, xj):
    sum = 0
    for xi in range(W.shape[0]):
        prod = getProd_except_j(i, j, messages, xi)
        psi = getPsi(xi, xj) 
        phi = getPhi(xi)
        total_prod = prod * psi * phi
        sum = sum + total_prod
    return sum

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


def updateMessage(i, j, xj, messages, new_message):
    #print("!!Message update!!")
    messages[xj][i][j] = new_message
    #messages[xj][j][i] = new_message #maybe problem??????????
    #print("New matrix: ")
    #print(messages)

def getNormalConst(messages, i):
    running_sum = 0
    rows, cols = A1.shape
    for xi in range(W1.shape[0]):
            running_sum += (getPhi(xi)*getProd(i, xi, messages))
    nc = running_sum
    return nc

def getBelief(i, messages):
    b_arr = np.zeros((W1.shape[0]))
    normal_const = getNormalConst(messages, i)
    for xi in range(W1.shape[0]):
        phi = getPhi(xi)
        #print("Phi   == ", phi)
        prod = getProd(i, xi, messages)
        #print("Prod  == ", prod)
        bi_xi = (phi * prod)/normal_const
        #print("Bi_Xi == ", bi_xi)
        #print()
        b_arr[xi] = bi_xi
    return b_arr

def getPairwiseBeliefs(messages, pw_beliefs):

    for xi in range(W1.shape[0]):
        for xj in range(W1.shape[0]):
            phi_i = getPhi(xi)
            phi_j = getPhi(xj)
            psi_ij = getPsi(xi, xj)
            prod_i = getProd_except_j(i, j, messages, xi)
            prod_j = getProd_except_j(j, i, messages, xj)


def getColoringAssignment(beliefs):
    asgnmt_list = []
    for belief in beliefs:
        x = np.argwhere(belief == np.amax(belief))
        x = x.flatten().tolist()
        if(len(x)>1):
            color = 0
        else:
            color = x[0]+1
        #x = np.argmax(belief)
        #print("For belief: ", belief, " x = ", x)
        asgnmt_list.append(color)
    return asgnmt_list

# def getPartitionFunc(beliefs):


def sum_product(A, W, its):
    r, c = A.shape
    messages_t = initializeMessages(A, W)
    messages_t_minus_1 = initializeMessages(A, W)
    beliefs_t = initializeBeliefs()
    # pw_beliefs_t = initializePwBeliefs()
    #print("Initial message matrix: \n", messages_t)
    #Do we need a copy of messages_t????????????????????????????? messages_t+1

    for t in range(its):
        #print("ITERATION ----------- ", t, " ----------------")
        for i in range(r):
            for j in range(c):
                if(A1[i][j]==1):
                    scaling_factor = getScalingFactor(W, messages_t, i, j)
                    for xj in range(W.shape[0]):
                        new_message_val = getNewMessage_SP(W, i, j, messages_t_minus_1, xj)
                        #scaling_factor = 1 #PROBLEM???????????????????????????????????????
                        #scaling_factor = getScalingFactor(xj, messages_t)
                        #print("SCALING FACTOR === ", scaling_factor)
                        scaled_message = new_message_val/scaling_factor
                        updateMessage(i, j, xj, messages_t, scaled_message) #REVIEW THIS
                        #xi = xj
                        #belief_xi = getBelief(xi)
            # beliefs__for_i = getBelief(i, messages_t)
            # beliefs_t[i] = beliefs__for_i
        for i in range(r):
            beliefs__for_i = getBelief(i, messages_t)
            beliefs_t[i] = beliefs__for_i

        #getPairwiseBeliefs()

        # print("!!@@!!@@!!@@--- Messages at t-1: ")
        # print(messages_t_minus_1)
        # print("!!@@!!@@!!@@--- Messages at t: ")
        # print(messages_t)
        messages_t_minus_1 = np.copy(messages_t)
    
    print("Belief matrix becomes (each row = 1 vertex, each col = color assignment): \n")
    print(beliefs_t)
    print()

    print("Hence coloring assignment for [v1, v2, ..., vn] becomes: ")
    assignment_list = getColoringAssignment(beliefs_t)
    print(assignment_list)

    #Z = getPartitionFunc(beliefs_t)

    # print("Message matrix becomes: ")
    # print(messages_t)
                        
def max_product(A, W, its):
    r, c = A.shape
    messages_t = initializeMessages(A, W)
    messages_t_minus_1 = initializeMessages(A, W)
    beliefs_t = initializeBeliefs()
    for t in range(its):
        for i in range(r):
            for j in range(c):
                if(A[i][j]==1):
                    #SCALING FACTOR
                    for xj in range(W.shape[0]):
                        new_message_val = getNewMessage_MP(W, i, j, messages_t_minus_1, xj)
                        #SCALED MESSAGE
                        updateMessage(i, j, xj, messages_t, new_message_val)
        for i in range(r):
            beliefs__for_i = getBelief(i, messages_t)
            beliefs_t[i] = beliefs__for_i

        messages_t_minus_1 = np.copy(messages_t)
    
    # print("Belief matrix becomes (each row = 1 vertex, each col = color assignment): \n")
    # print(beliefs_t)
    # print()

    
    assignment_list = getColoringAssignment(beliefs_t)
    #print(assignment_list)
    return assignment_list
    # print("Message matrix becomes: ")
    # print(messages_t)
    # print()


#-----------------------------------------------------------------------------------
#Graph 1: 1--2--3
#Weights: [1, 2]

global A1 
global W1 
# A1 = np.array([[0,1,0], [1,0,1], [0,1,0]])
# W1 = np.array([1,1])

# A1 = np.array([
#     [0, 1, 0, 0],
#     [1, 0, 1, 0],
#     [0, 1, 0, 1],
#     [0, 0, 1, 0]
# ])
# W1 = np.array([1, 2, 3])

A1 = np.array([
    [0, 1, 0, 0, 0],
    [1, 0, 1, 0, 0],
    [0, 1, 0, 1, 0],
    [0, 0, 1, 0, 1],
    [0, 0, 0, 1, 0]
])
W1 = np.array([1, 2, 3, 4, 5])
its1 = 30

user_input = "1"
while(user_input!=3):
    user_input = input("1) TYPE 1 for Sum-Product\n2) TYPE 2 for Max-Product\n3) TYPE 3 to Exit\n")
    
    if(user_input=="1"):
        sum_product(A1, W1, its1)

    elif(user_input=="2"):
        assignment_list = max_product(A1, W1, its1)
        print("------- MAX PRODUCT --------")
        print("Coloring assignment for [v1,v2,...,vn]: ", assignment_list)
        print("----------------------------\n")

    elif(user_input=="3"):
        print("Exiting...")
        break

    else:
        print("!!ERROR: Not a valid input")
    



#print("Your matrix: \n", A1)

