Index:

    Files
        sample.csv          |  dataset of known compounds up to 1868
        ex1TPR.png          |  example of periodic table representation (image)
        getTables.py        |  Script to find all unique Rs, and construct respective TPRs  
        getElems.py         |  Script to obtain a list of all unique elements present in dataset
        makeVecs.py         |  Script to convert compounds into vectors. Used in findR.py
        TPs.py              |  Datafile. This is where periodic tables can be defined (as dictionaries)
        scratch.py          |  File to store code produced but not longer used
        MPIgetTabs.py       |  Parallelized implementation of getTables.py


    Directories
        Data                |  Where to store all data I'm producing
        inf                 |  Directory to hold files relevant to the main report (LaTeX, etc)
        timeExperim         |  Contains files prepared for estimation of running time of code on whole DS
