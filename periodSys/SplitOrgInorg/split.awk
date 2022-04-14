#!/usr/bin/awk -f

# Split between organic and inorganic substances, based in following criteria:
# Criteria are:
# Is organic if:
# 	contain C + set of elements (see elemOrganic.txt)
# Is inorganic if:
#	contain C + set of elements (see elemIonrg.txt), in a regex below under `ino`
#	don't contain C

BEGIN{
FS="\t"
OFS=FS
ino="He:|Be:|Ne:|Al:|Ar:|Sc:|Ti:|V:|Cr:|Mn:|Fe:|Co:|Ni:|Cu:|Zn:|Ga:|Ge:|Kr:|Y:|Zr:|Nb:|Mo:|Tc:|Ru:|Rh:|Pd:|Ag:|Cd:|In:|Sn:|Sb:|Xe:|La:|Ce:|Pr:|Nd:|Pm:|Sm:|Eu:|Gd:|Tb:|Dy:|Ho:|Er:|Tm:|Yb:|Lu:|Hf:|Ta:|W:|Re:|Os:|Ir:|Pt:|Au:|Hg:|Tl:|Pb:|Bi:|Po:|At:|Rn:|Fr:|Ra:|Ac:|Th:|Pa:|U:|Np:|Pu:|Am:|Cm:|Bk:|Cf:|Es:|Fm:|Md:|No:|Lr:"
}

$2~"C:"{	# If contains Carbon

if(match($2,ino)){	# Also contains some element in `ino`
	inorg_sub[$0]
}
else{
	org_sub[$0]
}
next
}

{	# Else, don't contain carbon
	inorg_sub[$0]
}

END{
print "====== inorganic substances ======", length(inorg_sub)
for(i in inorg_sub)	print i

print "====== organic substances ======", length(org_sub)
for(i in org_sub)	print i
}
