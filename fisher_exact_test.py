import math
import scipy.stats as stats
import sys

USAGE_MESSAGE = "usage: python %prog [-f input_file] [-c class1,class2,class3,class4] [-h] [-tail]"

def getOptions():
	import optparse
	desc = 'Script for computing the Fisher\'s exact test'
	parser = optparse.OptionParser(USAGE_MESSAGE, description=desc)
	parser.add_option('-t', '--tail', action='store', type='string', dest='tail', default="", help='Use parameter\n\'less\' for lower p-value\n\'greater\' for higher p-value')
	parser.add_option('-c', '--classes', action='store', type='string', dest='classes', default="", help='Four classes comma seaparated')
        parser.add_option('-f', '--file', action='store', type='string', dest='inputFile', default="", help='Each line has a header and four classes space separated')
	parser.add_option('-i', '--index', action='store', type='string', dest='index', default="", help='Header index')
	parser.add_option('-s', '--suppress', action="store_true", dest='s', help='Suppress output')
	(options, args) = parser.parse_args()
	if options.tail: 
		tail = options.tail
	else:
		tail = None
	if options.classes: 
		classes = options.classes.split(',')
	else:
		classes = None
	return args, tail, classes, options.inputFile, options.s, options.index

def getTail(tail):
	if tail == "less":
                return "less"
        elif tail == "greater":
                return "greater"
        elif tail == None:
                return "two-sided"
        else:   
                return None

def computeFisher(filePath, sided=getTail(None)):
	genesDict = {}	
	with open(filePath) as f:
		for line in f:
			classes = line.replace(" ", "\t").replace("\n", "").split("\t")
        		gene = classes[0]
        		oddsratio, pvalue = stats.fisher_exact([[classes[1], classes[2]], [classes[3], classes[4]]], alternative=sided)
	        	genesDict[gene] = [oddsratio, pvalue, classes[1:]]
	return genesDict

def printGenesDict(genesDict):
	print "#Gene\tP-Value\tOddsRatio"
	for gene in genesDict:
		oddsratio = genesDict[gene][0]
		pvalue = genesDict[gene][1]
		classes = genesDict[gene][2]
		print gene+"\t"+str(pvalue)+"\t"+str(oddsratio)+"\t"+str(classes[0])+"\t"+str(classes[1])+"\t"+str(classes[2])+"\t"+str(classes[3])
    
if __name__ == '__main__':
	args, tail, classes, inputFile, suppress, index = getOptions()	
	tail = getTail(tail)
	if tail != None:
		if inputFile:
			filePath = inputFile #args[0]
			genesDict = computeFisher(filePath, tail)
			if genesDict != None:
				printGenesDict(genesDict)
		elif classes:
			oddsratio, pvalue = stats.fisher_exact([[classes[0], classes[1]], [classes[2], classes[3]]], alternative=tail)
			if suppress:
				if index:
					print index,
				print classes[0], classes[1], classes[2], classes[3], pvalue
			else:
				print "OR:", oddsratio
				print "p-Value:", pvalue
				if float(pvalue != 0.0):
					print "Score:", 0-math.log(pvalue, 10)
		else:
			print USAGE_MESSAGE
	else:
		print "Correct the tail value."
