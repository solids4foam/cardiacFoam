#!/bin/bash
#------------------------------------------------------------------------------
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     cardiacFoamChangeCopyright
#     adapted from solids4foamChangeCopyright, and
#     adapted from foamChangeCopyright in foam-extend-4.1
#
# Description
#     Updates the header of .H and .C files.
#
#     To run the script on multiple files, use the following command:
#
#     $ find . \( -name "*.C" -o -name "*.H" \) | xargs -n 1 "./cardiacFoamChangeCopyright.sh"
#
#------------------------------------------------------------------------------

# Check if the correct number of arguments are provided
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <file>"
    echo "Wrong number of arguments, expected 1 found $#"
    exit 1
fi

# Check if the file exists
if [[ ! -f "$1" ]]; then
    echo "File not found: $1"
    exit 1
fi

# Read original file
orig_file=$(< "$1")

perl -0777 -p -i -e '
s/[ \t]+$//mg;

my $license = <<'\''LICENSE'\'';
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
LICENSE

my $banner_start = "/*---------------------------------------------------------------------------*\\\n";
my $banner_end = "\\*---------------------------------------------------------------------------*/\n";
my $full_header = $banner_start . $license . "\n" . $banner_end . "\n";

if (
    s{\A(/\*[-]+\*\\\n)(.*?)(\n\\\*[-]+\*/\n?)}{
        my ($start, $body, $end) = ($1, $2, $3);

        $body =~ s/\A\n+//;
        $body =~ s/\n+\z//;

        if ($body =~ /(?:\A|\n)License\n/s) {
            $body =~ s{(?:\A|\n)License\n.*?(?=(?:\n[A-Z][A-Za-z0-9]*\n)|\z)}{
                my $prefix = substr($&, 0, 1) eq "\n" ? "\n" : "";
                $prefix . $license
            }se;
        } elsif (length $body) {
            $body = $license . "\n\n" . $body;
        } else {
            $body = $license;
        }

        $body =~ s/\n{3,}/\n\n/g;
        $start . $body . $end
    }se
) {
    # Existing header updated in place.
} else {
    $_ = $full_header . $_;
}
' $1

# Read modified file
mod_file=$(< "$1")

# Compare original and modified file
if [ "$orig_file" != "$mod_file" ]; then
    echo "File header does not match the expected header"
    echo "The header has been updated!"
    exit 1;
else
    echo "File header matches the expected header."
    exit 0;
fi
