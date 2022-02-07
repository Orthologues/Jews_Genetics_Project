#! /usr/bin/bash
find -name "*.main.txt" -mindepth 2 -maxdepth 2|while read file; do
    prefix=$(echo $file|sed -rE "s/\.main\.txt$//") &&
    pop_name=$(echo $prefix|cut -d $'/' -f 3|cut -d "_" -f 2) &&
    python parsing_fastGT_output.py --input $prefix --pop $pop_name \
    --output fastGT_summary --surr-pop-info newdata_180pops_info.csv
done
