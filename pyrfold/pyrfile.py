"""This module contains all of the scripts needed to create and read
files which store processed information that has been converted from
kinefold output
TODO
- consider removing sub files all together
- rework sub_file function to have a single dictionary or something just
easier to work with
- rework all of the submission structures to be less dictionary heavy
"""
########################################################################
############################ Modules ###################################
########################################################################
import os
import glob
import csv

########################################################################
############################   Read  ###################################
########################################################################
def sub_file(inputfile, justexperimentalconditions=False):
    """(str) -> {str:tuple}, {str:list}, {str:tuple}
    This should read the contents of the csv and present all of the information
    contained in two dictionaries

    >>> py_design_input_reader('test.csv')
        ({'TestName1': ('ACDGADDDDDD', 1, 2),
        'TestName2': ('ABSFDSADSFA', 1, 17),
        'TestName3': ('ABSFDSADSFADS', 1, 17)},
        {'TestName1': [['apple', 1, 3], ['zinc', 1, 2]],
        'TestName2': [['time', 1, 2]],
        'TestName3': [['time', 1, 2], ['pep', 1, 2]]},
        {'TestName1': (12,2),
        'TestName2': (113,3),
        'TestName3': (12,5)})
    """
    if justexperimentalconditions:
        # conditiondict = {}
        with open(inputfile, 'rU') as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None) #Skips the header
            for row in reader:
                if row[0]:
                    # conditiondict['polrate'] = float(row[4])
                    # conditiondict['dwelltime'] = float(row[5])
                    # conditiondict['fiveprime'] = int(row[2])
                    # conditiondict['threeprime'] = int(row[3])
                    # return conditiondict
                    return [row[4], row[5], row[2], row[3]]
    #initalizing variables:
    devicetosequence = {}
    devicetopart = {}
    devicetokinefoldparms = {}
    #Dictionary containing entanglements, cotrans, psudoknots
    devicetoexperimentalparms = {}
    devicetoforced = {}
    #Exctracting data
    with open(inputfile, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None) #Skips the header
        for row in reader:
            if row[0]:
                #defining the part name:(seq,lowbound,highbound)
                devicetosequence[row[1]] = (row[0].upper(),
                                            int(row[2]),
                                            int(row[3]))
                # polymerization, time after elongation (sec)
                devicetokinefoldparms[row[1]] = (float(row[4]), float(row[5]))
                #grabbing all of the part entries
                #They are stored in the end of the sheet
                #There are three things needed for device name,left,right
                # position
                devicetoexperimentalparms[row[1]] = (float(row[6]),
                                                float(row[7]), float(row[8]))
                #Collect forced Helix data
                forcedhelixes = []
                for index in range(9, 18, 3):
                    if row[index]:
                        forcedhelixes.append([int(row[index]),
                                              int(row[index + 1]),
                                              int(row[index + 2])])
                templist = []
                for i in range(18, len(row) - 1, 3):
                    if row[i]:
                        templist.append([row[i],
                                        int(row[i+1]),
                                        int(row[i+2]),
                devicetosequence[row[1]][0][int(row[i+1]) - 1: int(row[i+2])]])
                        #Last step adds the sequence of the part to the list
                        #1 is subtracked to account from the 1 to 0 index
                        #it is only subtracked from the low bound because of how
                        #python counts [n:c] calls
                devicetopart[row[1]] = templist
                devicetoforced[row[1]] = forcedhelixes
    return devicetosequence, devicetopart, devicetokinefoldparms, \
             devicetoexperimentalparms, devicetoforced

def sub_summary(filename):
    """(filename)->dict(partname:window)
    this should open the subsummary file and extract all the needed
    information for putting together the processed

    sub summary file stores start and stop of the parts within the
    shifted sequence as would be read with index 1
    IE needs to be shifed to work properly with python
    """
    devicetoparts = {}
    #open the .csv file
    with open(filename, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None) #Skips the header
        for row in reader:
            templist = []
            maxrange = len(row) - 1 #this is to avoid index errors
            for i in range(7, maxrange, 3):
                if row[i]: #This looks through all of the possible parts
                    templist.append([row[i], int(row[i+1]), int(row[i+2])])
            devicetoparts[row[0]] = templist
    return devicetoparts

def summary_file_to_dict(filename):
    """(filename)->dict(reference stucture name:[[sequence,freq],...])
    This will be useful for all comparisons of structure

    2013-12-11 10:42 WEV making changes to this to account for the new
    (..) notation that exists in these files


    """
    # Open the file
    outdictionary = {}
    with open(filename, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # Check to see if on name row or partlist row
            #name,
            #partname, sequence_n, frequencey_n, ...
            if (len(row) > 1) and row[1]: #should imply non name and non seq
                part, sequencelist = grep_from_partlist(row)
                sequence = next(reader)[0]
                outdictionary[part] = [sequence, sequencelist]
    return outdictionary

def summary_refpart_to_sequence(filename):
    """(filename)->dict(reference stucture name:[[sequence,freq],...])
    This will be useful for all comparisons of structure

    2013-12-11 10:42 WEV making changes to this to account for the new
    (..) notation that exists in these files
    """
    # Open the file
    outdictionary = {}
    with open(filename, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # Check to see if on name row or partlist row
            #name,
            #partname, sequence_n, frequencey_n, ...
            if (len(row) > 1) and row[1]: #should imply non name and non seq
                part, sequencelist = grep_from_partlist(row)
                outdictionary[part] = sequencelist #actually a [[],[],...]
    return outdictionary

def summary_exppart_to_sequence(filename):
    """(filename)->dict(concstruct name:{partname:[[sequence,freq],...]})
    This will be useful for all comparisons of structure
    READS IN the processed experimental summary
    """
    # Open the file
    outdictionary = {}
    tempdictionary = {}
    firstpass = 1
    with open(filename, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # Check to see if on name row or partlist row
            #name,
            #partname, sequence_n, frequencey_n, ...
            if (len(row) > 1) and row[1]: #should imply non name
                part, sequencelist = grep_from_partlist(row)
                tempdictionary[part] = sequencelist #actually a [[],[],...]
                #Let's just skip the sequence to make this workable for now
                next(reader, None)
            else: #if no data has been collected
                if firstpass:
                    constructname = row[0]
                    #print constructname + '!'
                    firstpass = 0
                else:
                    outdictionary[constructname] = tempdictionary
                    #print tempdictionary
                    tempdictionary = {}
                    constructname = row[0]
                    #print constructname
            #accounting for last entry
            outdictionary[constructname] = tempdictionary
    return outdictionary

def grep_from_partlist(sequencelist):
    """(list)->string, [[],[],...]
    the list of lists will contain [frequency, sequence]
    this is done for easy sorting

    >>> x = [['test', 'sequence1', 1, 'sequence2', 2, '', '', 'sequence3', 3]]
    >>> grep_from_partlist(x)
    >>> ('test', [['sequence1', 1], ['sequence2', 2], ['sequence3', 3]])
    """
    iterator = iter(sequencelist)
    #Grab first entery
    name = iterator.next()
    outlist = []
    for i, it in enumerate(iterator):
        #print it
        if it: #if not empty
            if not i % 2: #if first entry
                outlistcomponent = [float(it)]
            else: #If second entry
                outlistcomponent.append(it)
            if len(outlistcomponent) == 2: #only two entries required
                outlist.append(outlistcomponent)
    return name, outlist

def get_round_summaries_data(filename):
    """2014-01-10 WEV
    This is being modified to it's most basic form assuming that only
    filenames as well as frequences are stored in this file
    """
    outlist = []
    with open(filename, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None) #skips the ehader
        for row in reader:
            templist = [float(row[1]), row[0]]
            outlist.append(templist)
    return outlist

########################################################################
############################   Write  ##################################
########################################################################
def submission(devicetosequence, devicetopart, devicetokinefoldparms,
                    devicetoexperimentalparms, directorypath):
    """(dict, dict, dict) -> csv containing the summary of what
    was requested by the submission document
    """
    filepath = os.path.join(directorypath, 'sub_summary.csv')
    f = open(filepath, 'wb')
    #create header
    headers = ['name', 'windowed sequence', 'Polymerization Rate (nt/s)',
               'Folding time after elongation (s)', '1 renaturation 2 contrans', 'psudoknots 0 no 1 yes', 'entanglements 0 no 1 yes']
               #'sequence from experiment', 'size of sequence']
               #will try to incldue this information in a summary file
    for i in range(5):
        headers.append('part' + str(i+1))
        headers.append('part start')
        headers.append('part stop')
     #writing headers
    for name in headers:
        f.write(name + ',')
    f.write('\n')
    #shift the part windows to the window sequence
    for device, devicedata in devicetosequence.iteritems():
        #initializing bounds
        #devicetopart[device] should be a [[],[]]
        seq = devicedata[0]
        start = int(devicedata[1])
        stop = int(devicedata[2])
        seq = seq[start - 1: stop]
        #write name of device
        f.write(device + ',' + seq)
        f.write(',')
        f.write(str(devicetokinefoldparms[device][0]) + ',' +
                str(devicetokinefoldparms[device][1]) + ',')
        #Write experiment information
        for experimentalparms in devicetoexperimentalparms[device]:
            f.write(str(experimentalparms) + ',')

        for part in devicetopart[device]:
            #name, start shift, stop shift
            f.write(part[0] + ',')
            f.write(str(int(part[1]) - int(devicedata[1]) + 1) + ',')
            f.write(str(int(part[2]) - int(devicedata[1]) + 1) + ',')
            # shifted = [part[0], int(part[1]) - int(devicedata[1]) + 1,
            #            int(part[2]) - int(devicedata[1]) + 1]
        #     keytoadd.append(shifted)

        # windowdevicetosequence[device] = keytoadd
        f.write('\n')
    f.close()

def filled_in_form(filename, dicty):
    """()-> csv file to fill in
    This should simply create an empty submission file for submission
    dicty = {name: KineSubData class}
    """
    with open(filename + '.csv', 'wb') as f:
        headers = ['sequence', 'name', 'window start',
                   'window stop', 'Polymerization Rate (nt/s)',
                   'Folding time after elongation (s)',
                   '1 renaturation 2 contrans',
                   'psudoknots 0 no 1 yes',
                   'entanglements 0 no 1 yes']
                   #'sequence from experiment', 'size of sequence']
                   #will try to incldue this information in a summary file
        for i in range(5):
            headers.append('part' + str(i+1))
            headers.append('part start')
            headers.append('part stop')
       #writing headers
        for head in headers:
            f.write(head + ',')
        f.write('\n')
        for name in dicty:
            f.write(dicty[name].sequence + ',') #sequence
            f.write(name + ',') #name
            #windowStartStop[]
            f.write(str(dicty[name].windowstart) + ',' +
                                        str(dicty[name].windowstop) + ',')
            #polrate
            f.write(str(dicty[name].polrate) + ',')
            #folding time after elongation
            f.write(str(dicty[name].foldtimeafter) + ',')
            #Contrans or renaturation
            f.write(str(dicty[name].experimenttype) + ',')
            #psuedonots
            f.write(str(dicty[name].pseudoknots) + ',')
            #entanglments
            f.write(str(dicty[name].entanglements) + ',')
            partnamelist = dicty[name].partnamelist
            partstartstoplist = dicty[name].partstartstoplist
            for index, partname in enumerate(partnamelist):
                startstop = partstartstoplist[index]
                f.write(partname + ',' + str(startstop[0]) + ',' +
                                                    str(startstop[1]) +',')
            f.write('\n')

def blank_form(filename):
    """()-> csv file to fill in
    This should simply create an empty submission file for submission
    """
    f = open(filename + '.csv', 'wb')
    headers = ['sequence', 'name', 'window start',
                   'window stop', 'Polymerization Rate (nt/s)',
                   'Folding time after elongation (s)',
                   '1 renaturation 2 contrans',
                   'psudoknots 0 no 1 yes',
                   'entanglements 0 no 1 yes']
               #'sequence from experiment', 'size of sequence']
               #will try to incldue this information in a summary file
    forcedlist = ['forced start', 'forced stop', 'forced size']
    for i in range(3):
        headers.extend(forcedlist)
    for i in range(5):
        headers.append('part' + str(i+1))
        headers.append('part start')
        headers.append('part stop')
   #writing headers
    for name in headers:
        f.write(name + ',')
    f.write('\n')
    f.close()
