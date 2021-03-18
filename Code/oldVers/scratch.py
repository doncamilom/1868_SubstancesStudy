#! /usr/bin/env python3
# In this file I'll dump anything that isn't essential for execution (but may have use further on)

#########
#Implementation of first variant
#########

def getRs1(cmpnds,TP):
    # Find all unique Rs
    n = len(cmpnds[0])
    Rs = []
    ns = []
    
    for c in cmpnds:
        indx = np.nonzero(c)  #Get index of non-zero entries
        for i in indx[0]:     #Loop through these
            c_ = c.copy()
            ns.append(c_[i])  #Save subindex for this element. A unique R comes with a subind. e.g. R-Xn
            c_[i] = 0
            Rs.append(c_)     #Append the compound with a zeroed entry (R-X0)
    
    Rs = np.array(Rs)
    ns = np.array(ns).reshape(-1,1)
    
    Rs = np.concatenate([Rs,ns],axis=1)   # First columns represent an R ligand. last column represents subind. n
    Rs = np.unique(Rs,axis=0) #This array contains all unique Rs.
    print("Unique Rs: ", Rs.shape[0])
    
    # Now for each R(n), need to find all elements X such that compound R-Xn exists in dataset. 
    # Build a list of these for each R(n)
    
    R_list = []
    # Now generate the R representations on TP (TPR).
    j=0 #Counter for Rs with more than one appearence
    for R_ in Rs:
        n = R_[-1]  #Take subindex 
        R = R_[:-1] #The actual R
        curr_list = []  #List to hold elements X, where R-X exist, for this R
        
        for comp in cmpnds:  #Loop through all compounds
            elemInCmpnd = np.nonzero(comp)[0]
            for i in elemInCmpnd:  #loop through elements present in this compound
                tmp = np.zeros(comp.shape[0])        
                tmp[i] = n  
                if np.all((comp - tmp)==R) and comp[i] == n: #If removing element i from this compound gives the same R, then R-i exists 
                    curr_list.append(elemDict[i])  #Add this element to the list
        
        if len(curr_list)>1:  #Only consider Rs that exist in more than 1 compound
            np.save(f'TPR1/R{j}.npy', getTable(TP,curr_list)) 
            R_list.append(curr_list)
            j+=1
    
    return R_list  #This list contains lists (one for each R) of elements X such that R-X exist in dataset.
