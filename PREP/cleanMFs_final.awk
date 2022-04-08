#! /usr/bin/awk -f

# Produce the final tsv file containing usable MFs, with ID and date.
# This program removes any problematic (and unhandled) substances. Also removes repeated MFs coming from isomers.
# Modify this if you find ways around the special cases where MFs disagree.

BEGIN{
FS="\t"
OFS=FS
}

# Case: No differences between formulas
$3~0{
if($4){	# In case a MF exists indeed
	if($4 in y){
		if($2<y[$4]){
			y[$4]=$2
			id[$4]=$1
		}
	}
	else{
		y[$4]=$2
		id[$4]=$1
	}
next	# Next makes sure following block works as an 'else' statement of this block
}
}

# Other cases:	Use only stoichiometric substances (defined as substances without a point anywhere in LSF)
!($4~/[0-9]*\./){

s=gensub("<sub>|</sub>","","g",$6)		# Clean Hill formula: substitute C<sub>2</sub>... for C2...
if(s in y){
       if($2<y[s]){
       	y[s]=$2
       	id[s]=$1
       }
}
else{
       	y[s]=$2
       	id[s]=$1
}
}

END{	# Print only unique formulas
for(i in y)	print id[i], i, y[i]
}
