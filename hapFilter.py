import sys
import traceback
from gvcfLine import gvcfLine
from gvcfSample import gvcfSample
from gfilterMethods import gfilterMethods

if __name__ == "__main__":
    gvcf_loc = sys.argv[1]
    gvcf_out = sys.argv[2]

    outFile = open(gvcf_out,"w")

    CC_filter = gfilterMethods("CC_indMedDP.txt",12)
    GP_filter = gfilterMethods("GP_indMedDP.txt",14)

    with open(gvcf_loc) as gvcf:
        for raw_line in gvcf:
            try:
                fline = gvcfLine(raw_line)
                try:
                    if fline.isDataLine:
#                        if len(fline.alt) == 1 and len(fline.ref) == 1:
#                            outFile.write(raw_line)
#                        else:
#                            print raw_line
                        if CC_filter.siteFiltering(fline):
                            CC_samplePass = CC_filter.sampleFiltering(fline.CC_samples) 
                            GP_samplePass = GP_filter.sampleFiltering(fline.GP_samples)
    
                            if CC_samplePass ^ GP_samplePass:
                                 outFile.write(raw_line)
                            else:
                                print raw_line
                        else:
                            print raw_line
                    else:
                        outFile.write(raw_line) #Add header information

                except Exception, err:
                    print "failed data line/repr"
                    print fline.repr()
                    print (traceback.format_exc())
                    print raw_line

            except Exception, err:
               print "failed"
               print (traceback.format_exc())
               print raw_line

 
    gvcf.close()
    outFile.close()



