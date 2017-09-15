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
inputFolder="/share/data/tcga/coad/annovar_plus"
#inputFolder="/home/um99/tumors/COAD_broad"
parseAnnoutPath="/home/um99/scripts/parse_annout.py"
excludedPath="/home/um99/annovar_output/ContrastRunkPred/broad.excluded.txt"
tempPath="$path"
#effectMutationList="nonsynonymous_SNV"
effectMutationList="nonsynonymous_SNV,stopgain,stoploss,frameshift_insertion,frameshift_deletion"
geneColumn=6
normalColumn=10
primaryColumn=11
nameColumn=11 #it is the same of primaryColumn for coincidence
threshold="0.001"

#annotate each cancer file
samples=$(ls "$inputFolder/"*.gz | wc -l)
echo "Found $samples samples."
echo

i=0
for sample in "$inputFolder/"*.gz; do
	let i++
	tempFile="$tempPath/tempVCF.vcf"
	echo "Extraction sample $sample		($i - $samples)"
	zcat "$sample" > "$tempFile" 
	echo "Parsing normal genoytpe..."
	python "$parseAnnoutPath" "$tempFile" $normalColumn -t "$effectMutationList" -f $threshold > "$tempPath/$(basename $sample).NORMAL.vcf"; 
	sed -e "s/$/\t$(basename "$sample").NORMAL.vcf/" -i "$tempPath/$(basename $sample).NORMAL.vcf"
	echo "Parsing primary genotype..."
	python "$parseAnnoutPath" "$tempFile" $primaryColumn -t "$effectMutationList" -f $threshold > "$tempPath/$(basename $sample).PRIMARY.vcf"; 
	sed -e "s/$/\t$(basename "$sample").PRIMARY.vcf/" -i "$tempPath/$(basename $sample).PRIMARY.vcf"
	rm "$tempFile"; 
done

echo
#count the number of time in which a gene is mutated in Normal and in Primary
echo "Extracting the list of all mutated gene with the respective sample..."
cat "$tempPath/"*.vcf | awk '{print $'$geneColumn', $'$nameColumn'}' | sort -u > "$tempPath/mutatedGenesForSample.txt";
echo "Creating the list of all the genes mutated only in normal samples..."
grep "NORMAL.vcf" "$tempPath/mutatedGenesForSample.txt" | awk '{print $1}' | sort | uniq -c > "$tempPath/genesMutatedNormal.txt";
echo "Creating the list of all the genes mutated only in primary samples..."
grep "PRIMARY.vcf" "$tempPath/mutatedGenesForSample.txt" | awk '{print $1}' | sort | uniq -c > "$tempPath/genesMutatedPrimary.txt"; 

#create list of all mutated genes
echo "Creating the list of all mutated genes among the samples..."
awk '{print $1}' "$tempPath/mutatedGenesForSample.txt" | sort -u > "$tempPath/genesList.txt"
#rm "$tempPath/mutatedGenesForSample.txt"

echo
#print fisher test input
echo "Constructing the input file for the computing of the Exact Fisher's test..."

for gene in $(cat "$tempPath/genesList.txt"); do 
	tumorMutations=$(grep -E "(^|\s)$gene($|\s)" "$tempPath/genesMutatedPrimary.txt" | awk '{print $1} END {if(!NR) print 0}');
        tumorNoMutations=$(expr $samples - $tumorMutations);
	normalMutations=$(grep -E "(^|\s)$gene($|\s)" "$tempPath/genesMutatedNormal.txt" | awk '{print $1} END {if(!NR) print 0}'); 
	normalNoMutations=$(expr $samples - $normalMutations); 
	printf "$gene\t$tumorMutations\t$tumorNoMutations\t$normalMutations\t$normalNoMutations\n";
done > "$tempPath/fisherInput.txt"
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
