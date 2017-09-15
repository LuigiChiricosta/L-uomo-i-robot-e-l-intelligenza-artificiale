#!/bin/bash
if [ "$1" != "" ] && [ "$2" != "" ]; then
        path1="$1"
	path2="$2"
	
	format="0"
	var1=$(head -n1 "$path1" | awk '{print NF}'); 
	for i in $(seq 2 $var1); do
		format=$format",1.$i"; 
	done
	var2=$(head -n1 "$path2" | awk '{print NF}'); 
        for i in $(seq 2 $var2); do
                format=$format",2.$i";
        done
	join -a1 -a2 -e "0" -o "$format" <(sort -k1,1 "$path1") <(sort -k1,1 "$path2")
else
	echo "Use:"
	echo "$0 path1 path2"
fi
