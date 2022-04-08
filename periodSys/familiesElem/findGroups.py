#! /usr/bin/env python3

"""
Run CV-based algorithms for finding and refining PS groups,
out of similarity matrices and the optimized orderings found with GAs.
"""

# Add ../similarity to path, so we can import functions from simMat.py
import os,sys
path = os.path.join(os.path.dirname( __file__ ), '..', 'similarity')
sys.path.append(path)

from simMat import *
import pickle
import cv2
import multiprocessing as mp


def main():
    global NP,elemList, simMat_yr, min_yr, Indivs_yr, resultPath, dataPath
    # Load all similarity matrices
    dataPath = root+'Data/'
    resultPath = root+'Results/'
    elemList = getElemList(dataPath)
    simMat_yr = np.load(dataPath + 'simMat.npy')
    min_yr = 1800

    # Load the optimized orderings
    fh = open('./Results/optim_permut_yearly.gen', 'rb') 
    Indivs_yr = pickle.load(fh)

    with mp.Pool(processes=NP) as pool:
        results = [pool.apply_async(smoothForYear,
                                    args=(yr,10) )
                   for yr in range(1800,2022)]
        results_get = [r.get() for r in results]
    

    logger(results_get)
    filehandler = open(resultPath+'foundGroups.bin', 'wb')
    pickle.dump(results_get, filehandler)


def getElemList(dataPath):
    elemList = []
    with open("{}/ElementList.txt".format(dataPath),'r') as f:
        for line in f:
            elemList.append(line.strip())
    return elemList

def genToElem(gen):
    """func to convert a gen (list of positions) into sequence of reshuffled elements"""
    order = ['_' for i in range(103)]
    for i,idx in enumerate(gen):
        order[idx] = elemList[i]
    return order

def getGroups_CV(yr,ind,show=False,th=10,b=1,blur_sz=7,ups=5,max_grp=20):
    """Use computer vision algorithms to find blocks of 
    similar elements in the diagonal of similarity matrix"""
    # Obtain simetrized and reshuffled matrix using year=yr and ind-th ordering
    P = plot_simMat_yr(simMat_yr,yr,min_yr=min_yr,raw=False,scale=15,ordering=Indivs_yr[yr][ind],show=False,EL=elemList)

    # Set diagonal elements equal to avg of sorroundings, so edge is clearer
    idx = np.diag_indices(P.shape[0]-1)
    P[idx] = (P[idx[0],idx[1]+1] + P[idx[0],idx[1]-1])/2

    ## Start processing image
    
    # Convert matrix values into integers btwn 0-255
    img = np.rint(((P/np.max(P))*255)).astype(np.uint8)
    
    # Upsample: Increase num. pixels by duplicating each pixel locally
    US_FACT = ups # Upsample factor
    upsamp = cv2.resize(img,dsize=(img.shape[0]*US_FACT,img.shape[1]*US_FACT),
                        interpolation=cv2.INTER_NEAREST)  
    
    # Blur operation. Avg to fade some spurious edges
    blur = cv2.medianBlur(upsamp, blur_sz)
    
    # Add padding, so the image is framed and groups to the edges can be detected
    padd = cv2.copyMakeBorder(blur, 1,1,1,1, cv2.BORDER_CONSTANT, None, 0)

    # Detect edges with Canny algorithm
    canny = cv2.Canny(padd, th , th + b)
    
    # Find rectangles
    cnts = cv2.findContours(canny, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)[1]

    cp=upsamp.copy()
    final_cnts = []
    for c in cnts:  # Loop through found contours
        x,y,w,h = cv2.boundingRect(c)
        # Conditions: Square on the diagonal, below certain size, above certain size
        if x==y and w < max_grp*US_FACT and w > US_FACT:
            final_cnts.append(c)

    # Now let's extract the actual groups from the squares
    seq = genToElem(Indivs_yr[yr][ind])
    grps = {} # Do as a dict so repeated entries are not added
    
    covered = np.ones(P.shape[0])  # Array to store which positions have been covered
    for c in final_cnts:  # Loop through found rectangles
        x,y,w,h = cv2.boundingRect(c)
        st = seq[x//US_FACT:x//US_FACT+w//US_FACT]  # Add mins and maxs to include border of square too
        
        covered[x//US_FACT:x//US_FACT+w//US_FACT] = 0  # These positions are covered already
        # Add group only if it contains more than one element
        if len(set(st))>1:  
            grps[str(set(st))] = set(st)
            
    # Now form groups of consecutive elements that aren't covered by any square
    
    # Calculate the length of each of the "non-found" groups
    c,count=0,False
    lens=[]
    for d in np.diff(covered):
        if d==1:      c+=1; count=True
        if d==0 and count:  c+=1
        if d==-1: lens.append(c); c=0; count=False
            
    # Now extract all elements that weren't grouped, and extract groups based on the lengths above 
    p = np.array(seq[:P.shape[0]])[covered.astype(bool)]
    cum=0
    for l in lens:
        if l>1:
            st=p[cum:cum+l]
            grps[str(set(st))] = set(st)
        cum+=l
    
    return list(grps.values())

## Now the functions for average or smoothening groups. Noise reduction.
def avg1(grps):
    """Average out a list of groups. 
    Given a list of sets, return a unique set containing the elements 
    present in at least 50% of the input sets"""
    
    # Join them, then count number of times each element appears
    join = []
    for l in grps:        join += list(l)
    
    # Get list of uniq elems, then count number of appearances of each
    uniq = set(join)
    accept = len(grps)//2  # Acceptance threshold

    final_grp = []
    for e in uniq:
        if join.count(e) >= accept:  # If num appear > threshold, add to group
            final_grp.append(e)
            
    return set(final_grp)

def mainIter(InitState, NIters=1):
    """Implement the main iteration of algorithm. 
    Produce agreement sets among various CV parameters."""
    
    def M2(A1,A2):         
        A1,A2 = set(A1),set(A2)
        return len(A1.intersection(A2))/np.sqrt(len(A1)*len(A2))
    
    CurrState = InitState.copy()
    
    for Iter in range(NIters):
        mean_m2 = 0 # Keep track of mean similarity
        FinalState = []
            
        for g,Pg in enumerate(CurrState):  # For each pack in InitState
            Pg_smooth = {}          # Dict to store the updated version of InitState (dict so no repeats)
            for Ai in Pg:                  # For each group in pack Pg
                sim_grps = []  # Store most similar groups to Ai
                sim_grps_score = [] # Store similarity score (M2)

                # Loop over all packs, excluding Pg
                for i,pck in enumerate(CurrState): 
                    if i!=g:
                        # Best match found in this pack
                        most_sim,best_m2 = None,0

                        for Aj in pck: # Then for every group in such pack 
                            #Compare how similar Ai is to Aj
                            if len(Ai)==0 or len(Aj)==0: logger(Ai,Aj)
                            m2 = M2(Ai,Aj)
                            if m2>best_m2: most_sim = set(Aj); best_m2 = m2

                        if most_sim:  # If not None
                            sim_grps.append(most_sim)
                            sim_grps_score.append(best_m2)

                Avgd_Ai = avg1(sim_grps)
                if len(Avgd_Ai)>1:      Pg_smooth[str(sorted(Avgd_Ai))] = Avgd_Ai
                
            # Before appending to FinalState, clean repeated sets 
            FinalState.append(list(Pg_smooth.values()))  # Add smoothened Pg to Final state
            mean_m2 += np.mean(sim_grps_score)
            
        CurrState = FinalState.copy()  # Set current state to the newly formed state
        if Iter > 1:
            if mean_m2/len(FinalState) == meanM2: break
        meanM2 = mean_m2/len(FinalState)
        logger("mean M2 = {:.3f} after {} iters.".format(meanM2,Iter+1))
        
        
        
    def getUnique(FState):
        """Input: list of lists of sets  [ [{},{},...], [{},{},...], ...]
        Return a single list containing all unique sets"""
        Out = []
        for l in FState: Out += l # Concatenate all resulting groups
        # Filter for only unique
        Out_d = {}
        for l in Out:
            Out_d[str(sorted(l))] = l
        return list(Out_d.values())
        
    Return = getUnique(FinalState)
    TotalUniqInitState = getUnique(InitState)

    logger("Loop finished after {} iterations.".format(Iter+1))
    logger("\t ** InitState contained {} unique groups.".format(len(TotalUniqInitState)))
    logger("\t ** FinalState contains {}\n".format(len(Return)))
    return Return


def smoothForYear(yr,N_lt):
    """Run two instances of smoothing algorithm:
        1. For a range of parameters, run algo. to obtain list of groups attached to individual i
        2. Run algo. over all 20 lists of groups, one per individual.
        
    N_lt is just the number of different values of variable `lt` to try. The higher, the more reliable, but also more expensive.

    TO-DO: Find good combinations of parameters, so results are reproducible + CV is reliable"""
    
    InitState_year = []
    for ind in range(20):
        # Generate 3² different packs of groups, using different combinations of parameters
        InitState = []

        ups = 15
        ht = 40        

        # Make low values of lt (exponentially) more likely
        prob = np.exp(-np.arange(0,ht)/20)
        prob = prob/np.sum(prob) 
        for blur_sz_add in range(2,10,2): # let blur_sz = ups + 2, 4, 6, ... 
            for lt in np.random.choice(np.arange(0,ht),size=N_lt,p=prob):
                # ^^ A random sample of values, between 0 and ht-1 
                blur_sz = ups + blur_sz_add
                
                b = ht-lt
                th = lt
                
                gs = getGroups_CV(yr,ind,th=th,b=b,blur_sz=blur_sz,ups=ups)

                ### Check len of sets
                for i in gs:
                    if len(i)<2: logger(i)

                InitState.append(gs)

        # Now apply algorithm to smoothen sets
        logger("Smoothening ind = {}...".format(ind))
        s_ind = mainIter(InitState, NIters=10)
        InitState_year.append(s_ind)

    logger("Now processing results from all individuals")
    res = (yr,mainIter(InitState_year, NIters=10))

    fh = open(resultPath+'foundGroups/grp{}.bin'.format(yr), 'wb')
    pickle.dump(res, fh)
    
    logger("Done with yr {}".format(yr))
    return res



if __name__=="__main__":
    global root

    import argparse
    parser = argparse.ArgumentParser(description="Compute groups of elements using CV algorithms.")

    parser.add_argument("--NumProc","-N", required=True, type=int,
                        help="Number of processors to use for computations.")
    parser.add_argument("--inp_dir","-i", required=False, type=str, default="./",
                        help="Root directory for all computations")
    args = parser.parse_args()
    
    NP=args.NumProc
    root=args.inp_dir

    def logger(s,logfile=root+'log/families_cv.log'):
        with open(logfile, "a") as f:
            f.write(str(s)+"\n")

    from time import time
    t0=time()
    main()

    logger("Total time: {:.3f} s".format(time()-t0))
