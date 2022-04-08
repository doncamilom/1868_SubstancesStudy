#!/usr/bin/perl
use strict;
#use warnings;
use 5.010;

use Data::Dumper;

open (my $fh, '<:encoding(UTF-8)', $ARGV[0]);	# Load data file

while (my $line = <$fh>) {      # Read each line
    my @currline = split /\t/,$line;    # Split in fields sep. by \t. Store in currline
    print $currline[0], "\t";           # Print first field (ID)

    # Start splitting field 2 (composition)
    my @data = split /(?<!^)(?=[A-Z])/,$currline[1];    # Split into ['C','H2','Cl3']
    for my $i (@data) {                                 # For each element in formula 
        my @subdata = split /([A-z]+)([0-9]*)/,$i;      # Split each of these into ['Cl','3']
    
        if ($i =~ /[0-9]/){  print $subdata[1],":", $subdata[2], "-";}                  # If there is a number, print H:n
        else              {  print $subdata[1],":", 1, "-";          }                  # else, print C:1
    }
    print "\t", $currline[2];           # Finally print year.
}

