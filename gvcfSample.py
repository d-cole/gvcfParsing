import sys
import traceback

class gvcfSample:
    """
    Represents a sample from a vcf file
    If the sample is at an alt site only GT and DP info is present
    If the sample is missing missingInfo = True and no other info is loaded


    UNSURE ABOUT ALT READS AT REF SITES 
    """ 

    def __init__(self,i_sampleString,i_atAltSite,format):# -> None 
        """ 
        Instantiate a vcfSample object and load appropriate info 
        """ 
        self.format = format.split(":") 

        #Sample strings at alt sites have a different format than ref sites
        self.atAltSite = i_atAltSite
        self.sampleString = i_sampleString
        
        #Samples with missing info cannot be parsed
        self.missingInfo = False 
        self.GT = None

        self.altReads = 0
        self.refReads = 0
        self.otherReads = 0

        self.GQ = None
        self.PL = None
        self.DP = None

        self.PGT = None
        self.PID = None        

        #Parse string representation of sample
        try:
            self.__parseVals() 

        except:
            print i_sampleString
            traceback.print_exc()
            sys.exit()
                
    def __parseVals(self):# -> None
        """
        Parses the string representation of this sample
        """
        #Do not load info for samples with missing info
        if "./." in self.sampleString:
            self.missingInfo = True
            return
        
        #Split sample string into various info GT, AD, DP ...
        splitSample = self.sampleString.split(":")

        if self.atAltSite:
            #Load sample info from the format GT:AD:DP:GQ:PL
            self.GT = splitSample[0] 
            try:
                self.DP = int(splitSample[2])
            except:
                self.DP = 0
            self.altReads = int(splitSample[1].split(",")[1])
            self.refReads = int(splitSample[1].split(",")[0])
            self.otherReads = self.DP - (self.refReads + self.altReads)

            if len(splitSample) >= 5:
                #GT:AD:DP:PGT:PID or GT:AD:DP:GQ:PL GT:AD:DP:PGT:PID:PL
                if "PGT" in self.format:
                    self.PGT = splitSample[3]             
                if "PID" in self.format:
                    self.PID = splitSample[4] 
                if "GQ" in self.format:
                    self.GQ = float(splitSample[3])
                if "PL" in self.format:
                    self.PL = splitSample[4]
               
        else:#At ref site sample format is GT:DP

            self.GT = splitSample[0]
            self.DP = float(splitSample[1])

            #All ref sites seem to have no occurance of alt reads NOT SURE YET
            self.refReads = float(self.DP)
            self.altReads = 0
            self.otherReads = 0
    
    def getAltReads(self):# -> int
        """
        Returns number of alt reads for this sample
        """
        return self.altReads

    def getRefReads(self):# -> int
        """
        Return ref read counts from this sample
        """
        return self.refReads
    
    def getOtherReads(self):# -> int
        """ 
        Return nonmajor alt read counts from this sample
        """
        return self.otherReads

    def toString(self):# -> str
        """
        Return a string representation of this sample
        """
        return "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" \
            %(str(self.sampleString),\
            self.GT,\
            str(self.altReads),\
            str(self.refReads),\
            str(self.GQ),\
            self.PL,\
            str(self.DP),\
            str(self.missingInfo),\
            str(self.PGT),\
            str(self.PID),\
            str(self.missingInfo),\
            str(self.format),\
            str(self.sampleString),\
            str(self.atAltSite))
   
