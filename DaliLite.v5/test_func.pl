#!/usr/bin/perl

use strict;
use warnings;

# Function to retrieve the location of a program
sub which {
    my ($program) = @_;

    # Split the PATH environment variable into individual directories
    my @path_dirs = split /:/, $ENV{'PATH'};

    # Iterate through each directory to find the program
    foreach my $dir (@path_dirs) {
        my $file_path = "$dir/$program";

        # Check if the file exists and is executable
        if (-x $file_path) {
            return $file_path; # Return the full path
        }
    }

    return; # Return undefined if the program is not found
}

# Example usage:
my $program_name = 'mpirun'; # Change this to the program you want to find
my $program_location = which($program_name);

if ($program_location) {
    print "Location of $program_name: $program_location\n";
} else {
    print "Program '$program_name' not found in PATH\n";
}

