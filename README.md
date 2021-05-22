# Index 

/Preprocess/
- TPs.py            Definition of sample Periodic Table
- format.pl         Transform csv data into new format for next processing
- getAllRs.awk      Make all possible rewrittings of formulas.
- nonUniqRs.awk     Find all similarity relationships.
- runAll.sh         Run the three scripts above in sequence.
- dataToPy.py       Convert above data into numpy readable format.
- Data/             Dir. Contains results of various steps of preprocessing.
- datasets/         Dir. Contains original (toy) datasets in various formats
    - cleanDS.awk   Clean raw dataset from bad entries + fabricate test years.

/CoreCalcs/             
- TPs.py            Definition of sample Periodic Table
- playgrnd.ipynb    Notebook for exploration of in-depth aspects of data.
- Similarity.ipynb  NB with analysis of similarity matrices. Everything reduces to this one. Organize and publish
- simMat.py         Calculate similarity matrices yearly. Also contains support modules for above NB. 
- HND.py            Supporting functions for plotting HND.
- genetic.py        Implementation of GA for opt of 1D PT.
- Genetic/          Contains results of optimization with GA

/inf/                   All files regarding written report of activities and results.

/Data/                  Contains the real results, ran at Leipzig.
- Results/          Contains plots obtained from results on all data.
