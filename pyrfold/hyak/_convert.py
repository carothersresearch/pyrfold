"""
TODO write up docstrings for all of these objects
TODO write the technical procedure of what is actually happening in
this script
In general:
- This script looks for all of the intersections that exist between helixs
- If no intersections are found there are no pseudoknots
- If there are intersections the algrorithm removes the heix which is inovled
in the most intersections
- The interesections are recounted and the most intersection or rightmost
helix are removed
- this removed set is then given a second bracket type (i.e. pseudoknots)
- After this all of the new brackettype go through the same process
(loooking for pseudoknots within this pseudoknot)
- This continues until it's complete
"""

class KineRmnOutput:
    """Class to store .rnm file that comes from kinefold into useable
    data
    """
    def __init__(self, line):
        #This is a structure line
        self.helix = ''
        self.sequence = ''
        self.dotbracket = ''
        if "bases" in line:
            #This is cotrans line
            self.structure = line.split('     ')[0].strip()
            self.time = (float(line.split('     ')[1].
                split(' reached after ')[1].split(' ms')[0].strip()))
            self.freeenergy = (float(line.split(' | ')[1].
                split('kcal')[0].strip()))
            self.basenumber = (int(line.split(':')[len(line.split(':'))-1].
                split('over')[0].strip()))
            self.basetotal = (int(line.split(':')[len(line.split(':'))-1].
                split('over')[1].split('base')[0].strip()))
            self.type = 'cotrans'
        else:
            #This is a meltanneal thing
            self.structure = line.split('     ')[0].strip()
            self.time = float(line.split('     ')[1].split(' reached after ')[1].split(' ms')[0].strip())
            self.freeenergy = float(line.split(' | ')[1].split('kcal')[0].strip())
            self.basenumber = ''
            self.basetotal = ''
            self.type = 'meltanneal'

    def addHelix(self, helixline):
        self.helix = helixline.split('H')[0].strip()

    def generateDotBracket(self):
        #this uses the convert function to generate dotbracket
        self.sequence, self.dotbracket = \
                kine_to_dotbracket(self.structure, self.helix)

#Dictionary for all of the bracket types
BRACKETDICTIONARY = {0: ['[' ,']'], 2: ['<' ,'>'], 1: ['{', '}'],
                        3:["q", "p"], 4:["+", "-"], 5:['1', '0'], 6:['6','9'] }

def kine_to_dotbracket(totalsequence, totalhelix):
    helix = totalhelix.replace(' ', '')
    sequence = totalsequence.strip()
    #print helix
    #print sequence
    sequencelist, structurelist, helixpositionlistoflist = \
        reduce_kine_output(helix, sequence)
    helixindexlistoflists = generate_helix_index(helix, helixpositionlistoflist)
    listofknots, knotcountlist = knotintersection(helixindexlistoflists)
    listofknotknockouts = generate_knot_knockoutlist(listofknots, knotcountlist)
    dictofknots = {}
    bracketcount = 0
    for knot in listofknotknockouts:
        dictofknots[knot] = BRACKETDICTIONARY[bracketcount]
    while len(listofknotknockouts) > 0:
        #print dictofknots
        #print listofknotknockouts
        #raw_input()
        bracketcount += 1
        listofknots, knotcountlist = knotintersection(listofknotknockouts)
        listofknotknockouts = generate_knot_knockoutlist(listofknots, knotcountlist)
        for knot in listofknotknockouts:
            dictofknots[knot] = BRACKETDICTIONARY[bracketcount]

    #print listofknotknockouts
    sequence1, dotbracket = generate_dot_bracket(helixpositionlistoflist, helixindexlistoflists,sequencelist, dictofknots)
    #print totalsequence
    #print helixpositionlistoflist
    #print helixindexlistoflists
    # for i, base in enumerate(sequencelist):
    #     print str(i+1) + '   ' + base + '   ' + structurelist[i] + '   ' +  dotbracket[i]
    return generate_dot_bracket(helixpositionlistoflist, helixindexlistoflists
                                        ,sequencelist, dictofknots)

def generate_knot_knockoutlist(listofknots, knotcountlist):
    listofknotknockouts = []
    while sum(knotcountlist) > 0:
        # print 'knotcountlist'
        # print knotcountlist
        # print len(knotcountlist)
        # print 'listofknots'
        # print listofknots
        # print len(listofknots)
        # print listofknotknockouts
        # raw_input()
        numberofarcs = len(knotcountlist) - 1
        max_value = max(knotcountlist)
        max_index = knotcountlist.index(max_value)
        if max_value > 1:
            listofknotknockouts.append(listofknots[max_index])
            listofknots[max_index] = (0, 0)
            listofknots, knotcountlist = knotintersection(listofknots)
        elif knotcountlist[0] == knotcountlist[numberofarcs]:
            listofknotknockouts.append(listofknots[numberofarcs])
            listofknots[numberofarcs] = (0, 0)
            listofknots, knotcountlist = knotintersection(listofknots)
    return listofknotknockouts

def generate_dot_bracket(helixpositionlistoflist, helixindexlistoflists, sequencelist, dictofknots={}):
    previoushelixrightbound = 0
    dotbracket = ''
    for helixposition in helixpositionlistoflist:
        dotbracket = dotbracket + '.' * (helixposition[0] -
                                                    previoushelixrightbound )
        #look for helix index set that contains a value within the bounds of
        #helix position
        for helixpair in helixindexlistoflists:
            if helixposition[0] <= helixpair[0] <= helixposition[1]:
                if helixpair in dictofknots:
                    dotbracket = dotbracket + \
                    (dictofknots[helixpair][0]*(helixposition[1] -
                                                       helixposition[0] + 1))
                else:
                    #Have found the left bound (
                    dotbracket = dotbracket + ('('*(helixposition[1] -
                                                       helixposition[0] + 1))
            elif helixposition[0] <= helixpair[1] <= helixposition[1]:
                #do something with )
                if helixpair in dictofknots:
                    dotbracket = dotbracket + \
                    (dictofknots[helixpair][1]*(helixposition[1] -
                                                       helixposition[0] + 1))
                else:
                    dotbracket = dotbracket + (')'*(helixposition[1] -
                                                       helixposition[0] + 1))
        previoushelixrightbound = helixposition[1] +1
    if len(sequencelist) > len(dotbracket):
        dotbracket = dotbracket + '.'*(len(sequencelist) - len(dotbracket))
     #cycle through interactions to determine if brackets match up
    # for i, base in sequencelist:
    #     print ''.join(sequencelist)
    # print dotbracket
    # for i, x in enumerate(structurelist):
    #    print str(i), sequencelist[i], structurelist[i], helixlist[i], dotbracket[i]
    return ''.join(sequencelist), dotbracket

def reduce_kine_output(helix, sequence):
    sequencelist = []
    structurelist = []
    for i, x in enumerate(sequence):
        if x in ['G','U','A','C','g','u','a','c']:
            sequencelist.append(x)
        else:
            structurelist.append(x)
    structurelist.append('')
    #index all of the helix interactions
    #A helix interaction is some form of bracket
    helixpositionlistoflist = []
    if sequence[0] == '[':
        for i, struct in enumerate(structurelist):
            #go through all of the possible [ ] ^ or ' 'pdda
            if struct == '[':
                singlehelix = [i]
            elif struct == '^':
                singlehelix.append(i - 1)
                helixpositionlistoflist.append(singlehelix)
                singlehelix = [i]
            elif struct == ']':
                singlehelix.append(i - 1)
                helixpositionlistoflist.append(singlehelix)
    else:
        #go through all of the possible [ ] ^ or ' '
        for i, struct in enumerate(structurelist):
            if struct == '[':
                singlehelix = [i + 1]
            elif struct == '^':
                singlehelix.append(i)
                helixpositionlistoflist.append(singlehelix)
                singlehelix = [i + 1]
            elif struct == ']':
                singlehelix.append(i)
                helixpositionlistoflist.append(singlehelix)
    return sequencelist, structurelist, helixpositionlistoflist

def generate_helix_index(helixstring, helixpositionlistoflist):
    """ This goes through the helix strings and collects the
    index as well as the number associated with the interactions
    """
    # print helixstring
    helixlist = []
    p = 0
    for i, x in enumerate(helixstring):
        if (p) > 0:
            p += -1
            continue
        if x != '-':
            tempnumber = x
            while helixstring[i + p + 1] != '-':
                tempnumber = tempnumber + helixstring[i + p + 1]
                p += 1
                if helixstring[i + p] == "'":
                    break
            helixlist.append(tempnumber)
    helixindexlistoflists = []
    for position, index in enumerate(helixlist):
        if index:
            temphelixlist = []
            #append the first half of the interaction
            temphelixlist.append(helixpositionlistoflist[position][0] + 1)
            for position2, index2 in enumerate(helixlist):
                if (index + "'") == index2:
                    #store the last half of the interation
                    temphelixlist.append(
                        helixpositionlistoflist[position2][0] + 1)
                    #erase the prime version of the number
                    helixlist[position2] = ''
                    helixindexlistoflists.append(tuple(temphelixlist))
                    break
    return helixindexlistoflists
    #Reworking this concept
    #this will grab all of the helix numbers in order
    #this order will then be used to determine the proper numbering

def knotintersection(listlisthelix):
    #initialize the number of knots or 'intersections' that exist
    # print listlisthelix
    numberofknots = 0
    listofknots = []
    for index1 in listlisthelix:
        for index2 in listlisthelix:
            if ((index2[0] < index1[0] < index2[1]) and
               (index2[0] < index1[1] < index2[1])):
                #Not a pseudoknot because it's contained within another
                # print 'internal loop'
                # print index1
                # print index2
                continue
            elif ((index2[0] < index1[0] < index2[1]) and
                 (index1[1] > index2[1])):
                # print 'pseudoknot loop1'
                # print index1
                # print index2
                #pseudoknot
                numberofknots += 1
                listofknots.append(index1)
                break
            elif ((index2[0] < index1[1] < index2[1]) and
                 (index1[0] < index2[1])):
                # print 'pseudoknot loop2'
                # print index1
                # print index2
                #psuedoknot
                numberofknots += 1
                listofknots.append(index1)
                break
    # print "#knots - " + str(numberofknots)
    # print listofknots
    #Counting the number of intersections a single helix has
    #The idea her is to find the minimal number of deletions required for this
    #Process
    knotcountlist = []
    for knots in listofknots:
        # print knots
        #count number of intersections
        intersectioncount = 0
        for others in listofknots:
            nonknot = 0 #(((..((..<<.))..))).>>
            if knots == others:
                continue
            if knots[0] < others[0] < knots[1]:
                #this is an intersection
                intersectioncount += 1
                nonknot += 1
            if knots[0] < others[1] < knots[1]:
                #this is an intersection
                intersectioncount += 1
                nonknot += 1
            if nonknot == 2:
                intersectioncount += -2
        knotcountlist.append(intersectioncount)
    # print knotcountlist
    # print listofknots
    return listofknots, knotcountlist
