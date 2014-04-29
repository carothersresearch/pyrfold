class FoldingSubData(object):
    """This object is used to carry all of the information needed for hyak
    simulation submissions. This code should be extensible to other folding
    simulations as well.
    """
    def __init__(self, name, sequence, windowstart=0, windowstop=0,
        partstartstoplist=[], partnamelist=[], referencepart=None,
        forcedhelixes=[], polrate=30, foldtimeafter=1, experimenttype=2,
        pseudoknots=0, entanglements=0, numberofsimulations=10):
        #Name of the device
        self.name = name
        #Total sequence of the submitted sequence
        self.sequence = sequence
        #Window start and stop with 1 indexing
        self.windowstart = windowstart
        self.windowstop = windowstop
        #Collectin of parts and their corresponding names
        self.partstartstoplist = partstartstoplist
        self.partnamelist = partnamelist
        self.positionrefpart = referencepart
        #Defining fixed helix interactions for simulations
        #[leftstart, rightstart, length]
        self.forcedhelixes = forcedhelixes
        #Polymerization rate (nt/s) and foldtime after (s)
        self.polrate = polrate
        self.foldtimeafter = foldtimeafter
        #Experiment type 2 = cotrans, 1 = renaturation
        self.experimenttype = experimenttype
        #Controlling pseudoknots and entanglements for pseudoknots
        self.pseudoknots = pseudoknots
        self.entanglements = entanglements
        #The number of simulations to be done for a specific part
        self.numberofsimulations = numberofsimulations

    #A constructor of this class from a given entry in a sequence
    @classmethod
    def from_csv_sub_file_line(cls, listfromcsv):
        """This class method is used to initialize this data from a csv file
        for ease of transfer from python structures to submission files
        """
        name = listfromcsv[1]
        sequence = listfromcsv[0].upper()
        windowstart = int(listfromcsv[2])
        windowstop = int(listfromcsv[3])
        numberofsimulations = int(listfromcsv[4])
        partstartstoplist = []
        partnamelist = []
        if len(listfromcsv) > 20:
            for i in range(20, len(listfromcsv) - 1, 3):
                if listfromcsv[i]:
                    partnamelist.append(listfromcsv[i])
                    partstartstoplist.append([int(listfromcsv[i+1]),
                                              int(listfromcsv[i+2])])
        forcedhelixes = []
        if len(listfromcsv) > 11:
            for index in range(10, 19, 3):
                if listfromcsv[index]:
                    forcedhelixes.append([int(listfromcsv[index]),
                                          int(listfromcsv[index + 1]),
                                          int(listfromcsv[index + 2])])
        positionrefpart = listfromcsv[19]
        if not positionrefpart:
            positionrefpart = None
        polrate = float(listfromcsv[5])
        foldtimeafter = float(listfromcsv[6])
        experimenttype = int(listfromcsv[7])
        pseudoknots = int(listfromcsv[8])
        entanglements = int(listfromcsv[9])
        output = cls(name, sequence, windowstart, windowstop,
                    partstartstoplist, partnamelist, positionrefpart,
                    forcedhelixes, polrate, foldtimeafter, experimenttype,
                    pseudoknots, entanglements, numberofsimulations)
        return output

    def generate_csv_line(self):
        """Method to take the parameters that have been assigned and converting
        them into a csv line that can be written to a common csv sequence
        """
        outputlist = []
        outputlist.append(self.sequence)
        outputlist.append(self.name)
        outputlist.extend([self.windowstart, self.windowstop])
        outputlist.append(self.numberofsimulations)
        outputlist.append(self.polrate)
        outputlist.append(self.foldtimeafter)
        outputlist.append(self.experimenttype)
        outputlist.append(self.pseudoknots)
        outputlist.append(self.entanglements)
        for i in range(3):
            if self.forcedhelixes and len(self.forcedhelixes) - 1 >= i:
                outputlist.extend(self.forcedhelixes[i])
            else:
                outputlist.extend(['','',''])
        if self.positionrefpart:
            outputlist.append(self.positionrefpart)
        else:
            outputlist.append('')
        for counter, partname in enumerate(self.partnamelist):
            outputlist.append(partname)
            outputlist.extend(self.partstartstoplist[counter])
        return outputlist

    def kine_folding_data(self):
        "Kinefold requries total simulation time and pol rate in ms/nt"
        windowsize = self.windowstop - self.windowstart + 1
        if not self.polrate == 0:
            rate = float(1)/self.polrate*1000 #ms/nt
        else:
            rate = 0.
        requestedtime = rate*windowsize + self.foldtimeafter * 1000 #ms
        return rate, requestedtime

    def scale_parts_to_window(self):
        """This will rescale the part windows based on the window that was
        provided"""
        self.adjustedpartstartstop = {}
        start = self.windowstart
        stop = self.windowstop
        self.adjustedsequence = self.sequence[start-1:stop]
        for counter, partname in enumerate(self.partnamelist):
            tempstartstop = self.partstartstoplist[counter]
            tempstartstop[0] += -start + 1
            tempstartstop[1] += -start
            self.adjustedpartstartstop[partname] = tempstartstop

    def part_start_stop(self, partname):
        """This will return a the start and stop position of a part
        if it exists within the dictionary"""
        for counter, part in enumerate(self.partnamelist):
            if part == partname:
                return self.partstartstoplist[counter]
            else:
                print "Part " + partname + " not in partlist"








