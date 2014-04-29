
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
def folding_frequency(referencepartdict, experimentpartdict):
    """Compares experimentaldata to reference data to determine a realtive
    level of folding
    :param referencepartdict: A list of the referencsturcutres with freq for
        all parts
    :type referencepartdict: dict [part:[[freq, dotbracket], ...]
    :param experimentpartdict: A collection of all of the structures that
        appear for a given simulation run
    :type experimentpartdict: dict [part:['dotrbacket']:[dotrbacket:freq]]
    :return param: The collection of all the parts
    :return type: dict [part:frequency]
    """
    returndict = {}
    for part in experimentpartdict:
        returndict[part] = 0.
        if part not in referencepartdict:
            print "ERROR - part {} not in reference".format(part)
            break
        #Make the comparision
        for freqdotbracket in referencepartdict[part]:
            if freqdotbracket[1] in experimentpartdict[part]['dotbracket']:
                returndict[part] += experimentpartdict[part]['dotbracket'][freqdotbracket[1]]
    return returndict

def folding_frequency_one(refparttosequences,expparttosequences):
    """130829 WEV
    Custum script to generate the folding frequency data for the
    big insulating sequence run
    """
    # First loop through all devices:
    # concstruct name:{name stucture name:[[sequence,freq],...]}
    # build general solution dictionary for every experiment:
    exptofreqdict = {}
    for exp in expparttosequences:
        tempdict = {}
        for part in expparttosequences[exp]:
            freq = folding_part_in_any(refparttosequences[part],
                expparttosequences[exp][part], 0)
            tempdict[part] = freq
        exptofreqdict[exp] = tempdict
    return exptofreqdict

def folding_part_in_any(refstuctures, expstructures, tolerance):
    """(listoflist,listoflist,float)->float
    this function will compare all unique experimenal structures to the
    unique reference structures. If exp match any of the target folds
    it is seen as a success

    2013-12-11 16:20 listoflist = [[freq,dotbracket],[freq, dotrbacket]]
    This was changed to reflect this
    """
    output = 0
    for exppart in expstructures:
        for refpart in refstuctures:
            if refpart[1] in exppart[1]:
                #This adds the number of parts that were in the correct struct
                output += float(exppart[0])
                break
    return output

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
