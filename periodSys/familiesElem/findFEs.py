#! /usr/bin/env python3

"""
Run CV-based algorithms for finding and refining PS families of elements,
out of similarity matrices and the optimized orderings found with GAs.
"""

import sys
sys.path.append('similarity')
sys.path.append('Genetic1D')
from simMat import *
from genetic1D import *

import pickle
import cv2
import multiprocessing as mp

def main():
    global NP,elemList, simMat_yr, min_yr, Indivs_yr, resultPath, dataPath

    # Load all similarity matrices
    dataPath = rootPath + 'Data/'
    resultPath = rootPath + 'Results/'
    elemList = getElemList(dataPath)
    simMat_yr = np.load(dataPath + 'simMat.npy')
    min_yr = 1800

    # Load the optimized orderings
    fh = open(resultPath + 'optim_permut_yearly.gen', 'rb') 
    Indivs_yr = pickle.load(fh)


    with mp.Pool(processes=NP) as pool:
        results = [pool.apply_async(FEs_from_SM, args=(yr, 20))     # Use only 20 best permuts. for each year. 
                   for yr in range(1800,2022)]
        results_get = [r.get() for r in results]
    

    logger(results_get)
    filehandler = open(resultPath+'foundFEs.bin', 'wb')
    pickle.dump(results_get, filehandler)

######



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

def CVP(Y, ind, US_FACT, blur_sz, th, b, max_grp=20):
    """
    Computer Vision Pipeline developed for detection of families of elements in similarity matrices.
    Parameters:
        Y: Year to consider
        ind: index of the particular permutation to consider
        US_FACT: Upsampling factor
        blur_sz: Size of blurring kernel
        th: Hard Threshold (Canny's algorithm)
        b: Threshold With (Canny's algorithm)
    """

    # Get reshuffled matrix, using `ind`th permutation found for this year.
    P_reshuff = plot_simMat_yr(simMat_yr,Y,min_yr=min_yr,
                               raw=False,scale=15,
                               ordering=Indivs_yr[Y][ind],
                               show=False,EL=elemList)

    # Set diagonal elements equal to avg of sorroundings, so edge is clearer
    idx = np.diag_indices(P_reshuff.shape[0]-1)
    P_reshuff[idx] = (P_reshuff[idx[0],idx[1]+1] + P_reshuff[idx[0],idx[1]-1])/2

    # Convert matrix values into integers in range 0-255
    img = np.rint(((P_reshuff/np.max(P_reshuff))*255)).astype(np.uint8)

    # Upsample: Increase num. pixels by duplicating each pixel locally
    upsamp = cv2.resize(img,dsize=(img.shape[0]*US_FACT,img.shape[1]*US_FACT),
                        interpolation=cv2.INTER_NEAREST)

    # Blur operation. Avg to fade some spurious edges
    blur = cv2.medianBlur(upsamp, blur_sz)

    # Add padding, so the image is framed and groups to the edges can be detected
    padd = cv2.copyMakeBorder(blur, 1,1,1,1, cv2.BORDER_CONSTANT, None, 0)

    # Detect edges with Canny algorithm
    canny = cv2.Canny(padd, th , th + b)

    # Find rectangles
    cnts = cv2.findContours(canny, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)

    # ATTENTION: Depending on the version of OpenCV, contours are index 0 or 1.
    cnts = cnts[1]

    cp=upsamp.copy()
    final_cnts = []
    for c in cnts:  # Loop through found contours
        x,y,w,h = cv2.boundingRect(c)
        # Conditions: Square on the diagonal, below certain size, above certain size
        if x==y and w < max_grp*US_FACT and w > US_FACT:
            final_cnts.append(c)

    return final_cnts

def GetFEsCVP(Y, ind, final_cnts, US_FACT):
    """
    Post-process results from CVP. Convert raw info from squares, into FEs.
    """
    
    # Get similarity matrix for year Y, ordered with `ind`th permutation
    P = plot_simMat_yr(simMat_yr,Y,min_yr=min_yr,
                       raw=False,scale=15,
                       ordering=Indivs_yr[Y][ind],
                       show=False,EL=elemList)
    
    # Now let's extract the actual groups from the squares
    seq = genToElem(Indivs_yr[Y][ind])
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
    p = np.array(seq[:P.shape[0]])[covered.astype(bool)]

    c,count=0,False
    lens=[]
    for d in np.diff(covered):
        if d==1:      c+=1; count=True
        if d==0 and count:  c+=1
        if d==-1: lens.append(c); c=0; count=False

    cum=0
    for l in lens:
        if l>1:
            st=p[cum:cum+l]
            grps[str(set(st))] = set(st)
        cum+=l

    fams = list(grps.values())
    return fams

def Tanimoto(A1,A2):        
    A1,A2 = set(A1),set(A2)
    return len(A1.intersection(A2))/len(A1.union(A2))

def PopulMech(collect):
    """
    Popularity reduction of collection of FEs, into a single FE.
    Input: Collection of FEs
    Output: Single FE, whose elements exist in at least 50% of the FEs in input collection.
    """
    
    # Join them, then count number of times each element appears
    join = []
    for l in collect:   join += list(l)
    
    # Get list of unique elems, then count number of appearances of each
    uniq = set(join)
    accept = len(collect)*0.5  # Acceptance threshold

    final_fe = []
    for e in uniq:
        if join.count(e) >= accept:  # If num appear >= threshold, add to output FE
            final_fe.append(e)
            
    return set(final_fe)

def SNR(Collection, NIters, ID_run):
    """
    Implementation of Statistical Noise Reduction algorithm described.
        Takes as input a set of collections of FEs, and a max. number of iterations.
        Outputs a single collection of FEs.
    """
    CurrState = Collection.copy()
    
    for It in range(NIters):
        # Keep track of mean similarity between most similar FEs.
        sum_sim = 0 
        counter = 0  # Divisor for takin avg (sum_sim/counter)
        FinalState = []
        
        # Iter over collections.
        for i, Ci in enumerate(CurrState):
            Ci_smooth = {}
            
            # For each FE in Ci, get most similar FE in every Cj 
            for FE_x in Ci:
                sim_FEs = []
                sim_FEs_Tanim = []
                
                # Loop again over collections
                for j, Cj in enumerate(CurrState):

                    mostsim_FE, best_sim = None, 0
                    for FE_y in Cj:
                        sim = Tanimoto(FE_x, FE_y)
                        if sim > best_sim: 
                            mostsim_FE = set(FE_y)
                            best_sim = sim
            
                    if mostsim_FE:    # If best similarity > 0, append to list of FEs similar to FE_x
                        sim_FEs.append(mostsim_FE)
                        sim_FEs_Tanim.append(best_sim)
    
                # Get single representative FE
                condense = PopulMech(sim_FEs)
                if len(condense)>1:
                    Ci_smooth[str(sorted(condense))] = condense

                sum_sim += np.mean(sim_FEs_Tanim)
                counter += 1

            FinalState.append(list(Ci_smooth.values()))

            
        CurrState = FinalState.copy()  # Set current state to the newly formed state
        
        # Interrupt loop if meanSimIter converges to some value.
        if It > 1:
            if sum_sim/counter == meanSimIter: break
                
        meanSimIter = sum_sim/counter
        logger(f"ID {ID_run} Mean Similarity in {It+1}th Iteration = {meanSimIter:.3f}.")

        
    # All this is just diagnostic.
    def getUnique(FState):
        """
        Input: list of lists of sets  [ [{},{},...], [{},{},...], ...]
        Return a single list containing all unique sets
        """
        Out = []
        for l in FState: Out += l # Concatenate all resulting groups
        # Filter for only unique
        Out_d = {}
        for l in Out:
            Out_d[str(sorted(l))] = l
        return list(Out_d.values())
        
    Return = getUnique(FinalState)
    TotalUniqInitState = getUnique(Collection)

    logger(f"ID {ID_run} Loop finished after {It+1} iterations.")
    logger(f"ID {ID_run} \t ** InitState contained {len(TotalUniqInitState)} unique groups.")
    logger(f"ID {ID_run} \t ** FinalState contains {len(Return)}\n")

    return Return

def applyCVP_multipleParams(Y, ind, N_th=10):
    """
    Apply CVP using multiple combinations of parameters.
    Output a set of collections of FEs, each corresponding to some set of parameters.

    Input: 
        Y: year to study
        ind: opt permutation to use
        N_th: Number of different thresholds to sample for Canny.
    """

    Collections = []

    upsamp = 15

    hard_thresh = 40    # Same as th + b
    prob = np.exp(-np.arange(0,hard_thresh)/20)
    prob = prob/np.sum(prob) 
    th_list = np.random.choice(np.arange(0,hard_thresh), size=N_th, p=prob)

    for blur_w in range(2,10,2):
        blur_sz = upsamp + blur_w
        for th in th_list:

            gs = CVP(Y, ind, US_FACT=upsamp, blur_sz=blur_sz, th=th, b=hard_thresh-th)
            gs = GetFEsCVP(Y, ind, gs, US_FACT=upsamp)
            Collections.append(gs)

    return Collections


def FEs_from_SM(Y, Nind):
    """
    Apply the whole pipeline explained in Figure 1.
    Map from SM for year Y, all the way to a single collection of FEs.
    """

    # For this year, sort optimized permutations by cost.
    S = simMat_yr[Y - min_yr].copy()
    P = symmetrize(S)

    a = np.array([cost(P,ind) for ind in Indivs_yr[Y]])
    # Sort the list of Individuals for each year, based on cost
    order = np.argsort(a)

    collects = []

    for i in range(Nind):
        # Use ith best permutation
        InitState = applyCVP_multipleParams(Y, order[i])
        ci = SNR(InitState, NIters=15, ID_run=f"{Y}:{i}")
        collects.append(ci)


    logger("*** Now calculating final collection for year {}***\n".format(Y))

    res = ( Y, SNR(collects, NIters=15, ID_run=f"{Y}:FULL") )

    fh = open(resultPath+'foundFEs/FE_{}.bin'.format(Y), 'wb')
    pickle.dump(res, fh)
    
    logger("Done with year {}".format(Y))
    return res


if __name__=="__main__":
    global rootPath

    import argparse
    parser = argparse.ArgumentParser(description="Compute families of elements using CV algorithms.")

    parser.add_argument("--NumProc","-N", required=True, type=int,
                        help="Number of processors to use for computations.")
    parser.add_argument("--inp_dir","-i", required=False, type=str, default="./",
                        help="Root directory for all computations")
    args = parser.parse_args()
    
    NP=args.NumProc
    rootPath=args.inp_dir

    def logger(s,logfile=rootPath+'log/families_cv.log'):
        with open(logfile, "a") as f:
            f.write(str(s)+"\n")

    from time import time
    t0=time()
    main()

    logger("Total time: {:.3f} s".format(time()-t0))
