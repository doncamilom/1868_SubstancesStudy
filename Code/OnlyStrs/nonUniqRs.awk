#! /usr/bin/awk -f

BEGIN{OFS="\t"}
{
if ($1 in uniq_names){ # If R in uniq, then it's not unique.

    if (!($1 in nonuniq_names)){   # If this R is not in nonuniq, create an entry for it and add two first elems 

        nonuniq_names[$1];
        nonuniq[$1 "_elems"] = uniq[$1 "_elems"];
        nonuniq[$1 "_ID"]    = uniq[$1 "_ID"];
        nonuniq[$1 "_year"]  = uniq[$1 "_year"];
    }

    # Now add also this run's data
    nonuniq[$1 "_elems"] =  nonuniq[$1 "_elems"] ":" $2 ;
    nonuniq[$1 "_ID"]    =  nonuniq[$1 "_ID"] ":" $3 ;
    nonuniq[$1 "_year"]  =  nonuniq[$1 "_year"] ":" $4 ;

    }

else    # If R hasn't been observed, then add it to uniq
    uniq_names[$1];

    uniq[$1 "_elems"] =  $2 ;
    uniq[$1 "_ID"]    =  $3 ;
    uniq[$1 "_year"]  =  $4 ;

}


END{ 
    for (c in nonuniq_names) {
        print c,nonuniq[c "_elems"],nonuniq[c "_ID"],nonuniq[c "_year"]
    }
}
