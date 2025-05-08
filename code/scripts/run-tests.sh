#!/bin/bash

echo "========================================="
echo "TESTS"
echo "-----------------------------------------"

#!/bin/bash

cd tests || exit 1

# Define colors
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
END="\e[0m"

# Create an associative array to group files by steps-size
declare -A groups

# Collect files and group them
for file in *.bmp; do
    IFS='-' read -r prog threads steps size_ext <<< "$file"
    size="${size_ext%.bmp}"
    key="${steps}-${size}"
    groups["$key"]+="$file "
done

# Compare files within each group
ok=0
failed=0
for key in "${!groups[@]}"; do
    files=(${groups[$key]})
    echo -e "${YELLOW}Comparing files for steps-size: $key${END}"
    
    for ((i = 0; i < ${#files[@]}; i++)); do
        for ((j = i + 1; j < ${#files[@]}; j++)); do
            f1="${files[i]}"
            f2="${files[j]}"
            if cmp -s "$f1" "$f2"; then
                echo -e "  ${GREEN}[O] $f1 == $f2${END}"
                ok=$(($ok + 1))
            else
                echo -e "  ${RED}[X] $f1 != $f2${END}"
                failed=$(($failed + 1))
            fi
        done
    done
done


echo "========================================="
echo -e "[${GREEN}${ok} passed${END}/${RED}${failed} failed${END} tests]\n"