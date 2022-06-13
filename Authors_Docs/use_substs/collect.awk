BEGIN{
	FS="\t"
	OFS=FS
}


{
# count[yr][field]+=num
count[$1][$2]+=$3
}

END{
for(yr in count){
	for(field in count[yr]){
		print yr, field, count[yr]
	}
}

}
