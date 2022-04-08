#!/usr/bin/awk -f

# Group elements that are similar thorugh similarity relations, along with information from the corresponding compounds
# e.g. in H2O: OH-H , O-H2 , H2-O

BEGIN{ 
    FS="\t";
    OFS=FS;
}

{
e[$1]=e[$1] ":" $2
id[$1]=id[$1] ":" $3
yr[$1]=yr[$1] ":" $4
c[$1]+=1
}

END{
for(r in e){
	if(c[r]>1)	print	r, e[r], id[r], yr[r]
}
}

