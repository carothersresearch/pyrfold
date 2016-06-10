
"""This module exists to process the converted kinefold output into
useful metrics

Things TODO
- Change final structure to a specific time point
- Automatically calculate folding rates if desired
- Clean up sort and consolidatedotrbacket
    - I think I should be able to do this with a built in class
    I believe a set should do the trick

    I actually think a lot of this this could be solved using sets
    for comparisons etc...

"""
########################################################################
############################ Modules ###################################
########################################################################
def folding_frequency(referencepartdict, experimentpartdict, foldcutoff=0.0,
                      dominantstructure=False):
    """Compares experimentaldata to reference data to determine a realtive
    level of folding
    :param referencepartdict: a finalstructure dictionary like object that
        contains a list of parts with target structures
    :type referencepartdict: dict {'parts':{partname:Counter(dotbracket:freq)}}
    :param experimentpartdict: A collection of all of the structures that
        appear for a given simulation run
    :type experimentpartdict: dict {'parts':{partname:Counter(dotbracket:freq)}}
    :return param: The collection of all the parts
    :return type: dict [part:frequency]
    """
    returndict = {}
    #print referencepartdict
    for part in experimentpartdict['parts']:
        returndict[part] = 0.
        if part not in referencepartdict['parts']:
            print "ERROR - part {} not in reference".format(part)
            break
        #Make the comparision
        for expdotbracket in experimentpartdict['parts'][part]['dotbracket']:
            if dominantstructure and \
            (expdotbracket in
                referencepartdict['parts'][part]['dotbracket'].most_common(1)):
                returndict[part] += \
                 experimentpartdict['parts'][part]['dotbracket'][expdotbracket]
            elif expdotbracket in referencepartdict['parts'][part]['dotbracket']:
                #Check to see if the reference part frequency is greater than
                #The threshold
                if experimentpartdict['parts'][part]['dotbracket'][expdotbracket] > foldcutoff:
                    returndict[part] += \
                     experimentpartdict['parts'][part]['dotbracket'][expdotbracket]
    return returndict

def select_winners(resultsofpreviousround, foldfreqcuttoff):
    """2013-09-11 14:02
    results... = [[feqfold, nameofexperiment],...
    should order the array and then output the results
    """
    outlist = []
    resultsofpreviousround.sort(reverse = True)
    for exp in resultsofpreviousround:
        if exp[0] >= foldfreqcuttoff:
            outlist.append(exp[1])
        else:
            break
    return outlist
