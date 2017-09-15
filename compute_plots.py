import matplotlib.pyplot as plt
import numpy as NP
import math
import sys 

USAGE_MESSAGE = "usage: python %prog phdsnp_output_file -m -d [-h] [-b boxplot] [-f histogram] [-s saving]"
EFFECT="effect"
OCCURRENCE="occurrence"
PREDICTION="prediction"
SCORE="score"
POINT_CONSERVATION="point_conservation"
CONTEXT_CONSERVATION="context_conservation"
GENOTYPE="genotype"
FUNCTION="function"
columnDictionary = {EFFECT:8, OCCURRENCE:10, PREDICTION:12, SCORE:13, POINT_CONSERVATION:15, CONTEXT_CONSERVATION:16, GENOTYPE:17, FUNCTION:18}
PATHOGENIC = "Pathogenic"
BENIGN = "Benign"
GERMLINE = "germline"
SOMATIC = "somatic"
TSG = "TSG"
ONCOGENE = "oncogene"

def getColumnByName(column):
    try:
        return int(column)
    except:
        if column in columnDictionary:
            return columnDictionary[column]-1
        else:
            return None
            
def get_options():
        import optparse
        desc = 'Script for plotting distribution'
        parser = optparse.OptionParser(USAGE_MESSAGE, description=desc)
        parser.add_option('-m', "--mode", action='store', type='int', dest='m', help='mode (1: germline-somatic, 2: TSG-oncegene, 3: TSG+onogenes-not, 4:benign-pathogneic, 5:germline_benign-germline_pathogenic)')
        parser.add_option("-b", "--boxplot", action="store_true", dest='b', help='draw boxplot')
        parser.add_option("-f", "--histogram", action="store_true", dest='h', help='draw histogram')
        parser.add_option("-s", "--saving", action="store_true", dest='s', help='save plots')
        parser.add_option("-d", "--distribution", action="store", type='string', dest='d', help='column about distribution ("effect", "occurrence", "prediction", "score", "point_conservation", "context_conservation", "genotype", or any number)')
        parser.add_option('-t', "--snv-type", action='store', type='string', dest='t', help='snv type list comma separated ("synonymous_SNV, nonsynonymous_SNV, stopgain, stoploss, frameshift_insertion, frameshift_deletion")')
        parser.add_option('-x', "--x-axis", action='store', type='string', dest='x', help='x axis limits')
        parser.add_option('-y', "--y-axis", action='store', type='string', dest='y', help='y axis limits')
        parser.add_option("-l", "--log", action="store_true", dest='l', help='log scale')
        
        (options, args) = parser.parse_args()
        if options.d: 
            d = getColumnByName(options.d)
        else:
            d = None
        if options.t: 
            veff=options.t.split(',')
        else:
            veff=["nonsynonymous_SNV", "stopgain", "stoploss", "frameshift_insertion", "frameshift_deletion"]
        if options.x: 
            x=map(int, options.x.split(','))
        else:
            x=None
        if options.y: 
            y=map(int, options.y.split(','))
        else:
            y=None
        return args,options.m,options.b,options.h,options.s,d,veff,x,y,options.l

def computeKS(data1, data2):
    from scipy.stats import ks_2samp
    statistic, pValue = ks_2samp(data1, data2)
    return statistic, pValue

if __name__ == '__main__':
    argv = sys.argv;
    args,mode,boxplot,histogram,saving,distribution,veff,x,y,l=get_options()
    #inputPath = "C:\\Users\\Luigi\\Desktop\\Bioinformatics\\VB_shared_folder\\phdsnp.TCGA.output.genotype.vcf"
    outputPath = "C:\\Users\\Luigi\\Desktop\\Bioinformatics\\VB_shared_folder\\"
    error = False
    
    if len(argv) > 1:     
        if mode != None and distribution != None:
            effectID = getColumnByName(EFFECT)
            predictionID = getColumnByName(PREDICTION)
            genotypeID = getColumnByName(GENOTYPE)
            functionID = getColumnByName(FUNCTION)
            ditributionID = distribution#occurrenceID
            cutoff = False
            cutoffs = range(1, 100)
            effects = veff
            data1=[]
            data2=[]
            print "Extracting informations..."
            mutations = 0
            for data in open(argv[1]).readlines():
                mutations +=1
                dataSplit = data.replace(" ", "\t").split("\t")
                d = float(dataSplit[ditributionID].strip())
                e = dataSplit[effectID].strip()
                m = dataSplit[predictionID].strip()
                g = dataSplit[genotypeID].strip()
                f = ""
                if len(dataSplit) > functionID:
                    f = dataSplit[functionID].strip()
                #if cutoff == False or (cutoff == True and d in cutoffs):
                if e in effects:
                    if mode == 1: #the lists are either of germlines or somatics
                        if g == GERMLINE:
                            data1.append(d)
                        elif g == SOMATIC:
                            data2.append(d)
                    elif mode == 2: #the lists are either of TSGs or oncgenes
                        if f == TSG:
                            data1.append(d)
                        elif f == ONCOGENE:
                            data2.append(d)
                    elif mode == 3: #the lists are either TSGs,oncogenes,TSG/oncogenes or not
                        if TSG in f or ONCOGENE in f:
                            data1.append(d)
                        else:
                            data2.append(d)
                    elif mode == 4: #the lists are either of benigns or pathogenics
                        if m == BENIGN:
                            data1.append(d)
                        elif m == PATHOGENIC:
                            data2.append(d)
                    elif mode == 5: #the lists are both germlines either of benigns or pathogenics
                        if g == GERMLINE and m == BENIGN:
                            data1.append(d)
                        elif g == GERMLINE and m == PATHOGENIC:
                            data2.append(d)
                    else:
                        print "Mode error"
                        error = True
                        break
            
            if not error:                
                totalSize = len(data1)+len(data2)
                ratio1= len(data1)*100/totalSize
                ratio2= len(data2)*100/totalSize
                print "Size of input data: "+str(totalSize)+" (over "+str(mutations)+" mutations)\n- 1: "+str(len(data1))+" ("+str(ratio1)+"%)\n- 2: "+str(len(data2))+" ("+str(ratio2)+"%)"
                
                if len(data1) != 0 and len(data2) != 0:
                    statistic, pValue = computeKS(data1, data2)
                    print "Kolmogorov-Smirnov:\n- statistic:", statistic, "\n- pValue:", pValue
                else:
                    print "Kolmogorov-Smirnov non calcolabile."
                    
                if(boxplot):
                    print "Plotting boxplot..."
                    fig, ax = plt.subplots()
                    axes = plt.gca()
                    if y != None:
                        axes.set_ylim(y)
                    ax.set_title("Left => " + BENIGN + ": " + str(effects).replace("[", "").replace("]", "").replace("'", "") + "\n" + "Right => " + PATHOGENIC + ": " + str(effects).replace("[", "").replace("]", "").replace("'", ""))
                     
                    axPlot = [data1, data2]
                    #ax.boxplot(axPlot, 0, 'b+')
                    
                    bp = ax.boxplot(axPlot, 0, 'b+', patch_artist=True)
                    for box in bp['boxes']:
                        box.set(color='blue', linewidth=1)
                        box.set(facecolor = 'white')
                    for whisker in bp['whiskers']:
                        whisker.set(color='blue', linewidth=1)
                    for median in bp['medians']:
                        median.set(color='red', linewidth=1)
                    
                    effectString = ""
                    for effect in effects:
                        effectString += "."+effect
                    finalPath = argv[1] + effectString + ".boxplot"
                    
                    if(saving):
                        print "Saving in", finalPath
                        plt.savefig(finalPath + ".jpeg")
                        plt.savefig(finalPath + ".pdf")
                    print "Printing..."
                    plt.show()
                if(histogram):
                    print "Plotting histogram..."
                    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
                    ax1.set_title(str(effects).replace("[", "").replace("]", "").replace("'", ""))
                    ax2.set_title(str(effects).replace("[", "").replace("]", "").replace("'", ""))
                    
                    if x != None:
                        ax1.set_xlim(x[0], x[1])
                        ax2.set_xlim(x[0], x[1])
                    if y != None:
                        if l != None and y[0] == 0:
                            y[0] = 1
                        ax1.set_ylim(y[0], y[1])
                        ax2.set_ylim(y[0], y[1])
                    if l != None:# or y == None:
                        ax1.set_yscale('log', basey=10)
                        ax2.set_yscale('log', basey=10)
                    
                    if data1 != []:
                        bins1 = NP.array(range(int(math.floor(min(data1))), int(math.ceil(max(data1)))))
                        ax1.hist(NP.digitize(data1, bins1), bins1, alpha=0.7, linewidth=1, edgecolor = "black")
                    if data2 != []:
                        bins2 = NP.array(range(int(math.floor(min(data2))), int(math.ceil(max(data2)))))
                        ax2.hist(NP.digitize(data2, bins2), bins2, alpha=0.7, linewidth=1, edgecolor = "black")
                    
                    effectString = ""
                    for effect in effects:
                        effectString += "."+effect
                    finalPath = argv[1] + effectString + ".histogram.blablabla"
                    if(saving):
                        print "Saving in", finalPath
                        plt.savefig(finalPath + ".jpeg")
                        plt.savefig(finalPath + ".pdf")
                    print "Printing..."
                    plt.show()
        else:
            if distribution == None:
                print "Column for distribution not specified or not valid (-d)"
            if mode == None:
                print "Mode not specified (-m)"
    else:
            print USAGE_MESSAGE
