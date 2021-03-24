BEGIN{
    FS="\t"
}

{
    if (!match($2,/\./)){   # Remove non-stoichiometric compounds (any formula with a .)
        if ($2 in y){        
            if (y[$2] > $3) {   
                y[$2] = $3
                id[$2] = $1

                }
        }
        else { y[$2] = $3 ; id[$2] = $1 }
    }
}

END { for (i in y) print id[i], "\t", i, "\t", y[i] }

