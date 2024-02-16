package auxiliary_lib;

use strict;
use warnings;

# dbellieny: Function to retrieve the location of a program
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

1; # End of module
