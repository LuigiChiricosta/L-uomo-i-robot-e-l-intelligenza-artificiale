import scipy.stats as stats;
import sys

USAGE_MESSAGE = "usage: python %prog -f1 -c1 -f2 -c2"
KENDALL_TAU = "Kendall-Tau"
SPEARMAN = "Spearman"

def get_options():
        import optparse
        desc = 'Script for plotting distribution'
        parser = optparse.OptionParser(USAGE_MESSAGE, description=desc)
	parser.add_option("-k", "--kendall-tau", action="store_true", dest='k', help='Compute Kendall-Tau correlation')
	parser.add_option("-s", "--spearman", action="store_true", dest='s', help='Compute Spearman correlation')
        parser.add_option("-f", "--file1", action="store", type='string', dest='f1', help='First of file')
        parser.add_option('-g', "--file2", action='store', type='string', dest='f2', help='Second file')
        parser.add_option('-c', "--column1", action='store', type='int', dest='c1', help='Column that contains the values for the correlation in the first file')
        parser.add_option('-d', "--column2", action='store', type='int', dest='c2', help='Column that contains the values for the correlation in the second file')
        (options, args) = parser.parse_args()
        return args,options.f1,options.f2,options.c1,options.c2,options.k,options.s

def extractPositions(filePath, column):
	import operator
	dictionary = {}
	with open(filePath, "r") as correlationFile:
		for line in correlationFile:
			if line[0] != '#':
				lineSplit = line.replace(" ", "\t").split("\t")
				dictionary[lineSplit[0]] = float(lineSplit[column])
	return dictionary

def createList(dic1, dic2):
	list1 = []
	list2 = []
	for i in dic1:
		if i in dic2:
			list1.append(dic1[i])
			list2.append(dic2[i])
	return list1, list2

def selectCorrelation(value, name):
	if value != None:
		return name
	else:
		return None

def computeCorrelation(correlation, list1, list2):
	if correlation == KENDALL_TAU:
		return stats.kendalltau(list1, list2)
	elif correlation == SPEARMAN:
		return stats.spearmanr(list1, list2)
	else:
		return None

if __name__ == '__main__':
        argv = sys.argv;       
	correlationList = []
	args, file1, file2, column1, column2, k, s = get_options()
	column1 = column1-1
	column2 = column2-1
	correlationList.append(selectCorrelation(k, KENDALL_TAU))
	correlationList.append(selectCorrelation(s, SPEARMAN))
	if file1 != None and column1 != None and file2 != None and column2 != None:
		dic1 = extractPositions(file1, column1)
		dic2 = extractPositions(file2, column2)
		list1, list2 = createList(dic1, dic2)
		print "Length of file1:", len(dic1)
		print "Length of file2:", len(dic2)
		print "Length of final lists:", len(list1)
		for correlationName in correlationList:
			if correlationName != None:
				print correlationName, "correlation:"
				outcome = computeCorrelation(correlationName, list1, list2)
				print "Score:", outcome[0]
				print "p-Value", outcome[1]
       	else:
               	print USAGE_MESSAGE
