#!/usr/local/bin/perl
# Mike McQuade
# Median-dna.pl
# Given a set of DNA strings and an integer value
# of k, this program calculates the median string
# of all the DNA strings. If more than one k-mer is
# a median string, the program only returns one of them.

use strict;
use warnings;

# Initialize variables
my (@nums,$k,$medianString,$minHamming,@dna,@patterns);
my @alphabet = ('A','C','G','T');

# Open the file to read
open(my $fh,"<ba2b.txt") or die $!;

# Define variables with respective string and integers
@nums = split / /, <$fh>;
$k = $nums[0];
while (my $temp = <$fh>){
	chomp($temp);
	push @dna,$temp;
}
$minHamming = $k;

# Print out the median string
print Search();

# Close the file
close($fh) || die "Couldn't close file properly";

# Calculates all d-neighbors for the DNA string
sub dNeighbors {
	my (@dneighbors,@kmerArray);
	my $str = $_[0];
	my $d = $k;
	if ($d == 0) {return ($str);}
	for (my $iteration = 1; $iteration <= $d; $iteration++) {
		if ($iteration == 1) {push @kmerArray,$str}
		else {@kmerArray = @dneighbors;}
		
		foreach my $kmer (@kmerArray){
			for (my $i = 0; $i < length($kmer); $i++){
				foreach my $letter (@alphabet){
					my $distinct = $kmer;
					substr($distinct,$i,1) = $letter;
			
					# Only distinct strings are added to
					# the array. If the string is already
					# in the array, it does not get added
					# again.
					if (!grep(/^$distinct$/,@dneighbors)){
						push @dneighbors,$distinct;
					}
				}
			}
		}
	}
	return @dneighbors;
}

# Calculates all possible k-mers
sub allKmers {
	# Initialize variables
	my @words;
	my $currentPower = my $currentChar = 0;

	# Iterate through each place in the k-mer,
	# adding a value in each spot
	for (my $i = $k-1; $i >= 0; $i--) {
		$currentPower = scalar(@alphabet) ** $i;
		for (my $j = 0; $j < scalar(@alphabet) ** $k; $j++) {
			if ($i == $k-1) {push(@words,$alphabet[$currentChar]);}
			else {$words[$j] .= $alphabet[$currentChar];}

			# Decrement current power each time
			$currentPower--;

			# Increment current power each time
			if ($currentPower == 0){
				$currentChar = ($currentChar+1) % scalar(@alphabet);
				$currentPower = scalar(@alphabet) ** $i;
			}
		}
	}
	return @words;
}

# Searches the strings for the median string
sub Search {
	# Iterate through each possible k-mer
	foreach my $kmer (allKmers()) {			
		# Find the d-Neighbors of the given k-mer
		my @dArr = dNeighbors($kmer);

		# Check how many of the DNA strings the motif
		# or any of its d-Neighbors occur in.
		my $matches = my $distance = 0;
		foreach my $dnaStr (@dna){
			if (index($dnaStr,$kmer) != -1) {
				$matches++;
			} else {
				foreach my $slice (@dArr) {
					if (index($dnaStr,$slice) != -1) {
						$matches++;
						for (my $j = 0; $j < $k; $j++){
							# If there is a mismatch in the strings,
							# increment the Hamming distance by 1.
							if (substr($kmer,$j,1) ne substr($slice,$j,1)) {$distance++;}
						}
						last;
					}
				}
			}
			# If there exists a match in each string,
			# compare the calculated Hamming Distance
			# against the given lowest
			if ($matches == scalar(@dna)){
				# If the calculated distance is less
				# than or equal to the running minimum,
				# replace the minimum with the new distance
				# and replace the median string with the new k-mer
				if ($distance <= $minHamming){
					$minHamming = $distance;
					$medianString = $kmer;
				}
			}
		}
	}
	return $medianString;
}