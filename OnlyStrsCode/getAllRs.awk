#!/usr/bin/awk -f

BEGIN{ 
    FS="\t";
    OFS="\t";
}

{
    split($2,cmpnd,"-"); # Split compound into [elem1,n1,elem2:n2,...]

    for (elem in cmpnd) {       # Loop through elements in compound
        split(cmpnd[elem],tmp,":");   # Split into elem,n

        curr_elem=tmp[1]; # Current element
        curr_n   =tmp[2]; # Subindex of current element.

        split($2,refill,curr_elem ":" curr_n);           # refill: contains all before and all after current element

        counter=1;
        for (n_=curr_n-1; n_>=0; n_--) {  # Start decreasing in the range of this particular subindex. n, n-1, n-2, ...

            if (n_>0) {
                print refill[1] curr_elem ":" n_ refill[2] "X:" counter ,curr_elem, $1, $3;

            }
            else      {
                print refill[1] substr(refill[2],2) "X:" counter ,curr_elem, $1, $3;

            }
            counter++;
            
        }
    }
}
