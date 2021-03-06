#! /usr/bin/env python3

import multiprocessing as mp
import numpy as np
import sys

def main():
    global size
    N = int(sys.argv[1])
    size = int(sys.argv[2])

    print("Finding unique Rs...")
    Rlist = findRs(N)
    print("Done")


def getRepeated(Rs,i):
    print(f"Process {i}")
    return Rs
    
def findRs(N):
    Rs = np.random.random((N,200))
    n = 2
    ns = np.random.randint(1,n,size=(N,1))
    
    Rs = np.concatenate([Rs,ns],axis=1)

    # Create data chunks
    Rs_distrib_list = [Rs[Rs[:,-1]==i] for i in range(1,n+1)]

    with mp.Pool(processes=size) as pool:
        R_results = [pool.apply_async(getRepeated,
                                      args=(Rs_i,i))
                     for i,Rs_i in enumerate(Rs_distrib_list)]

        # blocks until all results are fetched
        R_results_get = [r.get() for r in R_results]   # A list of arrays

    # Merge filtered Rs by concatenating the resulting list
    Rs = np.concatenate(R_results_get,axis=0)

    return Rs

if __name__ == '__main__':
    main()


