#! /usr/bin/awk -f

# Compute similarity between a *given* substance and all over substances in dataset
# Takes as arguments: reference substance, and the subst_descriptor file

BEGIN{
	FS="\t";
	OFS="\t"


#line="27110009	4	-:110354:605301:605365:27165908:X-:102391:605283:742035:1730800:3587155:37229612:X-:103233:1730800:37229620:X-:110354:605301:605365:29370440:X	-38453444:1-56318574:1-56318575:1-42335773:1	-2014-2021-2021-2016"

#split(line,ref,FS)
#rts=ref[3]	# Set of reaction templates for reference substance
#
#split(rts,a,"-")
#for(i in a){mat=mat "|" a[i]}
#
#mat=substr(mat,3,length(mat))
#
}


## First, read substances
NR==FNR{i+=1;ref[i]=$0; next}

$1~ref[1]{
l1=$3
}

$1~ref[2]{
l2=$3
}

END{	#Calc sim
split(l1,a1,"-")
split(l2,a2,"-")

c=0
for(i in a1){
	for(j in a2){
		if(a1[i]==a2[j]){
			c+=1
			print a1[i],a2[j], c
		}
	}
}
print c

}


#$1~ge{
## Collect and concatenate information
#if(count[$2]<50000){	# Keep increasing string only if ammount of data is reasonable
#	templ[$2]=templ[$2] "-" $3
#	ids[$2]=ids[$2] "-" $4
#	dates[$2]=dates[$2] "-" $5
#	count[$2]+=1
#}
#else if(overcount[$2]==0){
#	# Activate flag and delete all information for this substance
#	overcount[$2]=1
#	templ[$2]=0
#	ids[$2]=0
#	dates[$2]=0
#	print "Ignoring substance with ID" , $2
#}
#}
#
#END{
#for(s in templ){
#	print s, count[s], templ[s], ids[s], dates[s]
#}
#}
