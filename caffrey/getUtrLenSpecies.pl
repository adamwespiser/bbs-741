#!/usr/bin/perl

use lib '.';
use ens;
use strict;
use warnings;

#my $sa_mouse=&getMouseSliceAdaptor();
#my $sa_human=&getHumanSliceAdaptor();
#&pOutUtrInfo($sa_mouse,"mouse");
#&pOutUtrInfo($sa_human,"human");

#print $ARGV[0]." is now being fetched....\n";

print "id,species,chr,transLength,threeUtrLength,fiveUtrLength\n";
&connect();
my @species = ("cat","chicken","chimpanzee","dog","fruitfly","human","mouse","platypus","rat","zebrafish");
map {&pOutUtrInfo(&getSpliceAdaptorSpecies($_),$_)} @species

