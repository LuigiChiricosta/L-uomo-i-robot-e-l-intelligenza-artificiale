#!/bin/bash
outputPath="/home/um99/annovar_output"
project="ContrastRunkProject"
path="$outputPath/$project"

if [[ -e $path ]] ; then
    i=1
    while [[ -e $path-$i ]] ; do
        let i++
    done
    path=$path-$i
fi

mkdir "$path"
echo "> Project folder path: "$(realpath "$path")
echo

fisherTestPath="/home/um99/scripts/fisherTest.py"
inputFolder="/share/data/1kgenomes/split"
parseAnnoutPath="/home/um99/scripts/parse_annout.py"
tempPath="$path"
tumorPath="/home/um99/annovar_output/ContrastRunkTCGA/ContrastRunkTCGANS"
geneColumn=6
normalColumn=10
primaryColumn=11
nameColumn=11 #it is the same of primaryColumn for coincidence
effectMutationList="nonsynonymous_SNV"
#effectMutationList="nonsynonymous_SNV|stopgain|stoploss|frameshift_insertion|frameshift_deletion"

#annotate each cancer file
samples=$(ls "$inputFolder/"*.vars | wc -l)
cancerSamples=$(ls "$tumorPath/"*PRIMARY.vcf | wc -l)
echo "Found $samples samples."
echo

echo "Adding file name to each file row..."
i=0
for sample in "$inputFolder/"*.vars; do
	let i++
	echo -ne "Working on sample $sample ($i/$samples)                   \r"
	cp "$sample" "$tempPath/" 
	sed -e "s/$/\t$(basename "$sample")/" -i "$tempPath/$(basename $sample)"
done

echo
#count the number of time in which a gene is mutated in Normal and in Primary
echo "Extracting the list of all mutated gene with the respective sample..."
cat "$tempPath/"*.vars | egrep -w "$effectMutationList" | awk '{print $'$geneColumn', $'$nameColumn'}' | sort -u > "$tempPath/mutatedGenesForSample.txt";
echo "Creating the list of all the genes mutated only in normal samples..."
awk '{print $1}' "$tempPath/mutatedGenesForSample.txt" | sort | uniq -c > "$tempPath/genesMutatedAmount.txt";

echo
#create list of all mutated genes
echo "Creating the list of all mutate genes among the samples..."
awk '{print $1}' "$tempPath/mutatedGenesForSample.txt" | sort -u > "$tempPath/genesList.txt"
#rm "$tempPath/mutatedGenesForSample.txt"

echo
#print fisher test input
echo "Constructing the input file for the computing of the Exact Fisher's test..."
i=0
genes=$(wc -l "$tempPath/genesList.txt" | awk '{print $1}')
for gene in $(cat "$tempPath/genesList.txt"); do
	let i++
	echo -ne "Working on gene $gene ($i/$genes)                       \r"
        tumorMutations=$(grep -E "(^|\s)$gene($|\s)" "$tumorPath/genesMutatedPrimary.txt" | awk '{print $1} END {if(!NR) print 0}');
        tumorNoMutations=$(expr $cancerSamples - $tumorMutations);
        normalMutations=$(grep -E "(^|\s)$gene($|\s)" "$tempPath/genesMutatedAmount.txt" | awk '{print $1} END {if(!NR) print 0}');
        normalNoMutations=$(expr $samples - $normalMutations);
        printf "$gene\t$tumorMutations\t$tumorNoMutations\t$normalMutations\t$normalNoMutations\n" >> "$tempPath/fisherInput.txt"
done #> "$tempPath/fisherInput.txt"
#rm "$tempPath/genesMutatedNormal.txt"
#rm "$tempPath/genesMutatedPrimary.txt"

echo
#compute fisher in one and two tail
echo "Computing the Exact Fisher's test (one tail)..."
python "$fisherTestPath" -f "$tempPath/fisherInput.txt" --tail "greater"> "$tempPath/fisherOutput1T.txt";
echo "Computing the Exact Fisher't test (two tail)..."
python "$fisherTestPath" -f "$tempPath/fisherInput.txt" > "$tempPath/fisherOutput2T.txt";
echo "Joining the results..."
join <(awk '{print $1, $2}' "$tempPath/fisherOutput1T.txt" | sort) <(sort "$tempPath/fisherOutput2T.txt") > "$tempPath/fisherOutput.txt"
#rm "$tempPath/fisherOutput1T.txt"
#rm "$tempPath/fisherOutput2T.txt"

echo
#final results
echo "Joining input and output file of the Exact Fisher's test..."
join <(sort "$tempPath/fisherInput.txt") <(sort "$tempPath/fisherOutput.txt") > "$tempPath/variantsOutcomes.txt"
#rm "$tempPath/fisherTestInput.txt"
#rm "$tempPath/fisherOutput.txt"

echo
echo "Finish."

echo
echo "See the output in the following path:"
echo "$tempPath/variantsOutcomes.txt"
echo
