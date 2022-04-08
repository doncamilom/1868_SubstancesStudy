#! /usr/bin/awk -f

# Extract Molecular Formulas (from substance data) and Year of publication (from reaction data)
# 

BEGIN{
# Here type in what fields you want to extract
a="IDE.LSF,IDE.MF,IDE.HSF,IDE.IMF,IDE.IHSF"

split(a,lookfor,",") # Array with fields to look for
FS="\t"
OFS=FS
}

# First load all dates from reaction data (extracted previously)
NR==FNR{
dates[$1]=$2
next
}

# Now get formulas
{
if(match($0,/<substance index/)){
	delete data  # Reinitialize array after a new substance field is found
	delete elc
}

if(match($0,/<IDE.XRN.+>/)){
	# Extract IDE XRN
	subsid=gensub(/.*<IDE.XRN.+>(.+)<\/IDE.XRN>/, "\\1","g",$0)
	}

# Collect data in each of the fields and append to array 'data'
for(i in lookfor){
	field=lookfor[i]
	if(match($0,field)){
		# Add content to data array
		data[field ":0"]+=1 # Increment vector's index
		rgx=".*<" field ">(.+)</" field ">"
		if(!match($0,/<YY.STR/))
			data[field ":" data[field ":0"]]=gensub(rgx, "\\1","g",$0)  # Write over this index
		}
}

if(match($0,"IDE.ELC")){
	# Add content to data array
	elc[field ":0"]+=1 # Increment vector's index
	rgx=".*<IDE.ELC>(.+)</IDE.ELC>"
	elc[field ":" elc[field ":0"]]=gensub(rgx, "\\1","g",$0)  # Write over this index
}


# Dump all the collected data into a single line
if(match($0,/<\/substance>/) && dates[subsid]){	# Output only if date exists
	line=""
	ref=data[lookfor[1] ":1"] # Set a reference to compare all against.
	discrep=0	# Assume no discrepancies between formulas so far

	for(i in lookfor){
		field=lookfor[i]
		for(j=1;j<=data[field ":0"];j++){	# Iterate over all formulas found
			line=line OFS data[field ":" j]	# concat
			if(ref!=data[field ":" j]) discrep=1	# If some field is found to be different from ref, report discrepancy	
		}
	}
	subl=""
	for(i=1;i<=elc[field ":0"];i++){	# Iterate over found elements of ELC, and concatenate
		subl=subl elc[field ":" i]
	}
	line=line OFS subl
	if(subl!=ref)	discrep=1

	if(discrep==0){
		ref=gensub("<sub>|</sub>","","g",ref)		# Clean formula: substitute C<sub>2... for C2...
		print subsid, dates[subsid], discrep, ref	# Print only one MF as they are all equal
	}
	else		print subsid, dates[subsid], discrep line	# Print all the information
}

# Done :)

}

