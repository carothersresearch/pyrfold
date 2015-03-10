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
import csv
from .foldingsub import FoldingSubData
import cPickle as pickle

########################################################################
############################   Read  ###################################
########################################################################


def load_pickled_sub_summary(picklefilepath):
    """The inputfile contains data that is stored as a pickled data structure
    this function simply loads that structure and generates additional
    information"""
    with open(picklefilepath, 'rb') as temppickle:
        outdict = pickle.load(temppickle)
    # Make sure the additional data is generated
    for device in outdict:
        outdict[device].scale_parts_to_window()
    return outdict


def load_pickled_reference_sub(picklefilepath):
    """he inputfile contains data that is stored as a pickled data structure
    this function simply loads that structure and generates additional
    information"""
    with open(picklefilepath, 'rb') as temppickle:
        outdict = pickle.load(temppickle)
    return outdict


def sub_file(inputfile, justexperimentalconditions=False):
    """ This will go through a submission file and create a dicitonary
    which contains FoldingSubData objecs which can be used to write
    all submission files
    """
    if justexperimentalconditions:
        # conditiondict = {}
        with open(inputfile, 'rU') as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None)  # Skips the header
            for row in reader:
                if row:
                    # conditiondict['polrate'] = float(row[4])
                    # conditiondict['dwelltime'] = float(row[5])
                    # conditiondict['fiveprime'] = int(row[2])
                    # conditiondict['threeprime'] = int(row[3])
                    # return conditiondict
                    return [row[4], row[5], row[2], row[3]]
    outdict = {}
    with open(inputfile, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # Skips the header
        for row in reader:
            if row:
                tempsubdata = FoldingSubData.from_csv_sub_file_line(row)
                outdict[tempsubdata.name] = tempsubdata
    return outdict


def sub_summary(filename):
    """(filename)->dict(partname:window)
    this should open the subsummary file and extract all the needed
    information for putting together the processed

    sub summary file stores start and stop of the parts within the
    shifted sequence as would be read with index 1
    IE needs to be shifed to work properly with python
    """
    devicetoparts = {}
    # open the .csv file
    with open(filename, 'rU') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # Skips the header
        for row in reader:
            templist = []
            maxrange = len(row) - 1  # this is to avoid index errors
            for i in range(7, maxrange, 3):
                if row[i]:  # This looks through all of the possible parts
                    templist.append([row[i], int(row[i+1]), int(row[i+2])])
            devicetoparts[row[0]] = templist
    return devicetoparts

########################################################################
############################   Write  ##################################
########################################################################


def filled_in_form(filename, devicenametosubobj):
    """()-> csv file to fill in
    This should simply create an empty submission file for submission
    dicty = {name: KineSubData class}

    TODO URGENT this is not working for the current submission system
    But it will work without forced data and parts
    """
    if '.csv' not in filename:
        filename += '.csv'
    with open(filename, 'wb') as f:
        writer = csv.writer(f)
        headers = ['sequence', 'name', 'window start',
                   'window stop', 'numberofsimulations',
                   'Polymerization Rate (nt/s)',
                   'Folding time after elongation (s)',
                   '1 renaturation 2 contrans',
                   'psudoknots 0 no 1 yes',
                   'entanglements 0 no 1 yes',
                   'helix minimum free energy (Leave Blank IF Unsure)']
        forcedlist = ['forced start', 'forced stop', 'forced size']
        for i in range(3):
            headers.extend(forcedlist)
        headers.append('posrefpart')
        headers.append()
        for i in range(5):
            headers.append('part' + str(i+1))
            headers.append('part start')
            headers.append('part stop')
        # writing headers
        writer.writerow(headers)
        if devicenametosubobj is None:
            pass
        else:
            for devicename in devicenametosubobj:
                linetowrite = devicenametosubobj[devicename].generate_csv_line()
                writer.writerow(linetowrite)


def blank_form(filename):
    filled_in_form(filename, devicenametosubobj=None)
