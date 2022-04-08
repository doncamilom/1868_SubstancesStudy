#!/usr/bin/awk -f


BEGIN{
FS="\t";
OFS=FS
}


{a[$2]+=1}

END{
print "Number of unique MFs", length(a); 

for(i in a){
	f[a[i]]+=1;
	if(a[i]>5000) print i,a[i]
};

print "#######";

for(i in f){
	print i, f[i]
}
}
