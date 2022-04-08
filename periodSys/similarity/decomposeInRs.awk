#!/usr/bin/awk -f

# Decompose Molecular Formulas into 'binary representations' that allow an explicit display of the role of elements in compounds
# e.g. in H2O: OH-H , O-H2 , H2-O

BEGIN{ 
    FS="\t";
    OFS=FS;
    MAXN=30
}

{
split($2,cmpnd,"-"); # Split compound into [elem1:n1,elem2:n2,...]

for (elem in cmpnd) {       # Loop through elements in compound
	split(cmpnd[elem],tmp,":");   # Split into elem,n
	
	curr_elem=tmp[1]; # Current element
	curr_n   =tmp[2]; # Subindex of current element.
	
	split($2,refill,curr_elem ":" curr_n);           # refill: contains all before and all after current element
	
	counter=1;
	for (n_=curr_n-1; n_>=0; n_--) {  # Start decreasing in the range of this particular subindex. n, n-1, n-2, ...
	
		# Generate line: entry of the form "C:6-H:3-N:3-O:3-X:2"
		if (n_>0){
			line=refill[1] curr_elem ":" n_ refill[2] "X:" counter;
		}
		else{
			line=refill[1] substr(refill[2],2) "X:" counter;
		}
		
		# Print line with all relevant information.
		if(counter<MAXN)	ind="eq " counter
		else			ind="gt " MAXN
		lines_vec[ind][line OFS curr_elem OFS $1 OFS $3]
		counter++;
		}
	}
}

END{
for(i in lines_vec){
	print "=============="
	for(l in lines_vec[i])	print l

}
}
