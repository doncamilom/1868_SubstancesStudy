#!/usr/bin/gawk -f

BEGIN{
	FS="\t"
	OFS=FS
}

$2{
y=$3
split($2, mf, "-")
for(cn in mf){
	split(mf[cn],cnlist,":")
	count[cnlist[1]][y] += 1
}
}

END{
# Print header: years
head="elem"
for(y=1800; y<=2022; y++)
	head=head OFS y
print head

for(e in count){
line=e
for(y=1800; y<=2022; y++){
	line=line OFS count[e][y]
}
print line

}
}

