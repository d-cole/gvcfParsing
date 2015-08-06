import sys
from gvcfLine import gvcfLine
from gvcfSample import gvcfSample
from gfilterMethods import gfilterMethods
"""
gvcfToCSV.py
Creates a .csv file with the variant information
"""

HETZ = "0/1"

COLUMNS = "sample,numRef,numAlt,GQ"

SAMPLES_CC = ["CC3-3_B","CC3-3_C","CC3-3_D","CC3-3_E","CC3-3_F",\
"CC3-3_H","CC3-3_I","CC3-3_J","CC3-3_K","CC3-3_L","CC3-3_N","CC3-3_O"]

SAMPLES_GP = ["GP2-3_B", "GP2-3_C", "GP2-3_D", "GP2-3_E", "GP2-3_F", "GP2-3_G",\
 "GP2-3_H", "GP2-3_I", "GP2-3_J", "GP2-3_K", "GP2-3_L", "GP2-3_M", "GP2-3_N", "GP2-3_O"]

def getMutantString(sampleList,namesList):
    csv_string = ""
    for i in range(0,len(sampleList)):
        if sampleList[i].GT  == HETZ:
            csv_string = namesList[i] + "," \
                + str(sampleList[i].refReads) + "," \
                + str(sampleList[i].altReads) + "," \
                + str(sampleList[i].GQ)
            break
    return csv_string


if __name__ == "__main__":
    file_name = sys.argv[1]
    outFile = open(file_name[0:-5] + ".csv",'w')
    outFile.write(COLUMNS + "\n")

    CC_filter = gfilterMethods("CC_indMedDP.txt",12) 
    GP_filter = gfilterMethods("GP_indMedDP.txt",14)

    with open(file_name) as f:
        for raw_line in f:
            line = gvcfLine(raw_line)
            if line.isDataLine:
                CC_string = getMutantString(line.CC_samples,SAMPLES_CC)
                GP_string = getMutantString(line.GP_samples,SAMPLES_GP)
                
                if CC_string != "" and GP_string != "":
                    CC_result = CC_filter.sampleFiltering(line.CC_samples)
                    GP_result = GP_filter.sampleFiltering(line.GP_samples)
                    
                    if CC_result and not GP_result:
                        outFile.write(CC_string + "\n")
                    if GP_result and not CC_result:
                        outFile.write(GP_string + "\n")
 
                if CC_string != "" and GP_string == "":
                    outFile.write(CC_string + "\n")  

                if GP_string != "" and CC_string == "":
                    outFile.write(GP_string + "\n") 

    f.close()
    outFile.close()     

