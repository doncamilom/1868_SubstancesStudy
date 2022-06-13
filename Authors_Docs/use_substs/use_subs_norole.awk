#!/bin/gawk -f

# Using this file, we can count how many times a substance has been used in a given year.
# ../../../Andres/reaction_template/Data/reaction_sets.tsv
BEGIN{
	FS="\t"
	OFS=FS

	SIDs="385737,1098229,969135,471223,1209228,605461,385801"
	split(SIDs,SIDl,",")	# List of substances to check
}

{
for(sid in SIDl){
	SID=SIDl[sid]
	if(match($3,SID)){
		count[SID][$4]+=1
	}
}
}

END{
for(sid in SIDl){
	SID=SIDl[sid]
	print "Substance ID = ", SID
	for(yr in count[SID]){
		print yr, count[SID][yr]
	}

}
}
