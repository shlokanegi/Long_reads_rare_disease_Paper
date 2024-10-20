import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import csv
import math
import random 

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

def getDataFileList(path="./site_beds"):
    #Get the list of bed files (modkit_comb methylbed files previously filtered to relevant sites)
    file_list=[]
    for root,dirname,filenames in os.walk(path):
        for file in sorted(filenames):
            if file.endswith("bed"):
                file_path= os.path.join(root,file)
                if not file.startswith("~"):
                     file_list.append(file_path)
    return file_list

def loadData(file_list=getDataFileList(path="./site_beds"), sig_file="data/CSS1.cpgs_sort2_p01.txt"):
    #loads all the methylation data from the filtered bed files
    #Generates the histogram
    #compares to the CSS1 episignature
    methPercent = {}
    siteSum = {}
    siteCount = {}
    for file_path in sorted(file_list):
        filename = os.path.basename(file_path)
        if filename.endswith(".bed"):
            sample = filename.split(".")[0]
            with open(file_path, "r", encoding="utf8") as sample_file:
                tsv_reader = csv.reader(sample_file, delimiter="\t")
                for row in tsv_reader:
                    (chrom,start,end,baseCode,score,strand,starp, endp, color, Ncov,pcMod, nmod, ncan,noth,ndel,nfail,ndiff,nnocall) = row
                    #Pulls out row from bed file and stores its percent methylation in dictionary
                    site = chrom+":"+end
                    methPercent[site] = methPercent.get(site,{})
                    methPercent[site][sample] = float(pcMod)
                    siteSum[site] = siteSum.get(site,0) + float(pcMod)
                    siteCount[site] = siteCount.get(site,0)+1
        
    for site in siteCount:
        #calculate the average methylation at the site for all samples
        siteSum[site] = siteSum[site] / siteCount[site]
        for sample in methPercent[site]:
            #Calculate difference from the mean methylation for each sample
            methPercent[site][sample] = methPercent[site][sample] - siteSum[site]
     
    
    with open(sig_file, "r", encoding="utf8") as ref_file:
        #Load the Episignature file (lifted over to hg38)
        tsv_reader = csv.reader(ref_file, delimiter="\t")
        for row in tsv_reader:
            (site, difMeth, pval) = row
            methPercent[site]["CSS1_signature"] = float(difMeth)*100
    
    scores = {}
    
    numRand=1000000
    randScores = [0 for i in range(numRand)]
    
    
    for site in methPercent:
        refVal = methPercent[site]["CSS1_signature"]
        
        for sample in methPercent[site]:
            #Calculate the fraction of sites where the signs of the CSS signature methylation difference (e.g. hyper vs hypomethylation) match the sign of the difference from mean methylation for each sample
            scores[sample] = round(scores.get(sample,0)+ refVal/math.sqrt(refVal*refVal)*methPercent[site][sample]/math.sqrt(methPercent[site][sample]*methPercent[site][sample])*.5/len(methPercent) + .5/len(methPercent),5)
            
        sample_list = [x for x in methPercent[site].keys() if x is not "CSS1_signature"]
        for x in range(numRand):
            #Calculate the fraction of sites where the signs of the CSS signature methylation difference (e.g. hyper vs hypomethylation) match the sign of the difference from mean methylation for each 1,000,000 permuted sample
            rv = methPercent[site][random.choice(sample_list)]
            randScores[x] = round(randScores[x] + refVal/math.sqrt(refVal*refVal)*rv/math.sqrt(rv*rv)*.5/len(methPercent) + .5/len(methPercent),5)

    randScores.sort()
   
    
    methPercent["score"] = {}
    methPercent["pv"] = {}

     
    for sample in scores:
        #Calculate the p-value for each non-random sample, e.g. the number of random samples as concordant or more concordant with the CSS1 episignature
        if sample is not "CSS1_signature":
            numGreater = len([x for x in randScores if round(x,2) >= round(scores[sample],2) ])
            pv = round(numGreater / (1+len(randScores)) * len(scores),2)
            if pv > 1:
                pv = 1
            methPercent["score"][sample]=scores[sample]
            methPercent["pv"][sample]=round(pv,2)
            print (sample + "\t" +str(scores[sample])+ "\t" + str(numGreater) + "\t" + str(pv) ) 
  
    
    
    methPercent["score"]["CSS1_signature"]=1
    
    dataframe = pd.DataFrame.from_dict(methPercent, orient='index')
    
    #Set up labels for the ouput heatmap
    for sample in scores:
        if sample is not "CSS1_signature":
           dataframe.rename(columns={sample: sample+" ["+str( round(methPercent["score"][sample],2) ) +"; p="+str( methPercent["pv"][sample] )+"]"}, inplace=True)
    
    
    dataframeSort = dataframe.sort_values(by="CSS1_signature", ascending=False).sort_values(by="score",ascending=True, axis=1)
    
    
    #add a blank column before the episignature column
    dataframeSort.insert(len(dataframeSort.columns) - 1,'',np.NaN)
    dataframeSort.drop("score",inplace=True)
    dataframeSort.drop("pv",inplace=True)


    #plot the heatmap    
    plt.figure(figsize=(28.8, 19.2))
    sns.heatmap(dataframeSort, cmap="RdYlGn", annot=True, fmt=".1f",yticklabels=dataframeSort.index,annot_kws={"size": 6})
    plt.title('Methylation Difference Heatmap', fontsize=25, fontweight='bold')
    plt.xlabel('Sample ID', fontsize=21, fontweight='bold', labelpad=20)
    plt.ylabel('CpG location', fontsize=21, fontweight='bold',labelpad=17)
    plt.subplots_adjust(top=0.92, right=1.05, left=0.13, bottom=0.16)

    heatmap_output_file = os.path.join("./", 'heat_map.svg')
    plt.savefig(heatmap_output_file, dpi=300)

    #save heatmap values to csv file
    output_file = os.path.join("./", 'heat_map.csv') #specify
    dataframeSort.to_csv(output_file, index=False)
    
    #plot the histogram
    plt.clf()
    plt.figure(figsize=(10, 6))
    sns.histplot(pd.DataFrame(randScores, columns=['scores']), binrange=[0.15,.95], stat='probability',bins=85,legend=False)
    plt.xlabel('Methylation Concordance with CSS1',fontsize=15)
    plt.ylabel('Frequency',fontsize=15)
    plt.title('Distribution of Permuted Samples',fontsize=18, fontweight='bold')

    plt.axvline(randScores[numRand-1], color='g', linestyle='--', linewidth=1, label='max')
    plt.axvline(randScores[0], color='g', linestyle='--', linewidth=1, label='min')

    
    plt.legend(loc='upper right',facecolor='w', edgecolor='black',framealpha=1)


    plt.annotate('PMGRC 146-146-0', xy=(scores["PMGRC-146-146-0"]-.002, 0.003), xycoords='data',
                     xytext=(-73, 25), textcoords='offset points',
                     arrowprops=dict(arrowstyle="->", color='r', lw=1.5),
                     bbox=dict(boxstyle="round,pad=0.3", edgecolor='r', facecolor='w'))

    last = -1
    n=0
    #Add X to histogram showing distribution of samples
    for sample in sorted(scores, key=scores.get):
        if sample is not 'CSS1_signature':
            score = scores[sample]
            n = n+.01/6
            if score > last:
                last = score
                n=0
            plt.annotate('x', xy=(score, n), xycoords='data',
                     xytext=(0, 1), textcoords='offset points', color='b')



    #Saves SVG of the histogram
    plt.savefig('score_histogram.svg', dpi=250)

    return dataframeSort



if __name__=="__main__":
    df = loadData()
