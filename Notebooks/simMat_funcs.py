#! /usr/bin/env python3

from Analysis import *

def main():
    global elemList
    
    dataPath = "./Data/"

    # Load element list
    elemList = []
    with open(f"{dataPath}/ElementList.txt",'r') as f:
        for line in f:
            elemList.append(line.strip())
    print(f"Number of elements: {len(elemList)}")
    sim_mat = calc_simMats_yearly(dataPath)
    np.save(dataPath+'history_simMat.npy',sim_mat)



def getSimilarities_yr(Tyr,element,mask,year):
    """Get array of the similarities between elements and the given element
    based on 'replaceability' in chemical formulas, for the CS existent at a given year."""
    # Filter out by year
    T = (Tyr <= year)&(Tyr>0)
    
    # Now start actually calculating simMats
    X,Y = TP[element] # Get coords of elem in PT
    # similarity: number of times other elements appear in same table
    simil = T[T[:,X,Y]].sum(axis=0)*1.
    
    # Apply mask: convert any "off table" entry into nan
    with np.errstate(invalid='ignore',divide='ignore'):
        simil /= mask
    # Extract similarity data for every element
    data_arr = {e:simil[TP[e][0],TP[e][1]] for e in TP.keys() if e in elemList}

    data_arr = pd.Series(data_arr).rename(element)
    return data_arr

def simMat(Tyr,mask,year):
    mat = pd.DataFrame(columns=[e for e in TP.keys() if e in elemList])
    for e in TP.keys():
        if e in elemList: # Calculate in order, but subjected to what's in elemList
            mat[e] = getSimilarities_yr(Tyr,e,mask,year)
    return mat

def calc_simMats_yearly(dataPath = "../OnlyStrsCode/Data/"):
    """Calculate all simMats yearly.
    First load all the data that was produced during preprocessing.
    Takes as input: dataPath, and number of elements in dataset"""

    # Load matches_pickle.bin
    Match_file = bz2.BZ2File(dataPath+'matches_pickle.bin', 'r')
    Matches = pickle.load(Match_file)
    miny,maxy = getMinMax(Matches)
    print(miny,maxy)

    Tyr = np.zeros((len(Matches),7,32))      # This gets the arrays filled with years
    for i,match in enumerate(Matches):
        getTable(TP, match[0], match[1], Tyr[i] )  # Match[1] to fill with years


    # mask: make 0 all entries where no element exists.
    mask = (Tyr>0).sum(axis=0)>0

    # Generate data yearly
    year_span = maxy-miny+1
    numElems = len(elemList)
    simMats_yr = np.zeros((year_span,numElems,numElems))
    
    for i,yr in enumerate(range(miny,maxy+1)):
        simMats_yr[i] = simMat(Tyr,mask,yr)

    return simMats_yr


if __name__ == '__main__':
    main()
