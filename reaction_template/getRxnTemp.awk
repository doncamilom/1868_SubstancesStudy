#! /usr/bin/awk -f

# This script extracts all possible reaction templates out of a list of sets of reactives.
# Takes as input the output from 'getRCTN.awk'


BEGIN{
	FS="\t";
	OFS="\t"}

{
#27084664:1	3	:1730800:17431054:761261	1981

n=split($3,rctns,":");
rxnid=$1;  # First elem of array is RXN ID


if(n>2){
# Loop thorugh reactants, and for each one create a sorted array with placeholder at the end
for(i=2;i<=n;i++){
	# Ignoring ith entry of rctns array

	delete rct_copy;
	c=1;
	for(j=2;j<=n;j++){
		if(j!=i){
			rct_copy[c]=rctns[j];  # We want to exclude ith entry, to extract reaction template it belongs to 
			#print  rctns[j];
			c++;
		}
	}

	asort(rct_copy) # Sort left reactants

	# Produce output string, containing description of reaction template
	st="";
	for(j=1;j<n;j++){
			st=st ":" rct_copy[j];  
	}
	st=st "X";

	# For the current substance in current reaction set, print substr, substance, rx_tmp, id, date
	# substr is last 4 digits of substance_id, for grouping later
	printf "%04i\t", substr(rctns[i],length(rctns[i])-3,length(rctns[i])) 
	print rctns[i], st, rxnid, $4 

	}
}

}


