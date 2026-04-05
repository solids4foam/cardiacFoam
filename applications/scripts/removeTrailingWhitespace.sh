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
#     removeTrailingWhitespace
#
# Description
#     Removes trailing spaces and tabs from tracked text files.
#
# Usage
#     ./applications/scripts/removeTrailingWhitespace.sh
#     ./applications/scripts/removeTrailingWhitespace.sh file1 file2 ...
#
#------------------------------------------------------------------------------

set -euo pipefail

if [[ $# -gt 0 ]]; then
    files=("$@")
else
    files=()
    while IFS= read -r -d '' file
    do
        files+=("$file")
    done < <(git grep -I -l -z '[[:blank:]]\+$' || true)
fi

if [[ ${#files[@]} -eq 0 ]]; then
    echo "No trailing whitespace found."
    exit 0
fi

echo "Removing trailing whitespace from:"
for file in "${files[@]}"; do
    echo "  $file"
done

for file in "${files[@]}"; do
    perl -0pi -e 's/[ \t]+$//mg' "$file"
done

echo "Done."
