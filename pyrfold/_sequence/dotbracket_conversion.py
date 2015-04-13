def reverse_complement_bracket_notation(bracket_notation):
    bracket_dict = {'(':')',
                    ')':'('}
    outstring = ''
    for base in bracket_notation:
        if base in bracket_dict:
            outstring += bracket_dict[base]
        else:
            outstring += '.'
    return outstring[::-1]


def reverse_numbering(helix_dictionary, size_of_dotbracket):
    outdict = {}
    for helix in helix_dictionary:
        outdict[helix] = []
        for start, stop in helix_dictionary[helix]:
            outdict[helix].append([size_of_dotbracket-start, size_of_dotbracket-stop])
    return outdict


from copy import copy


def dot_bracket_notation(bracket_notation):

    lastSquarredPairedPos = []
    lastRoundedPairedPos = []
    lastCurlyPairedPos = []

    basePairs = []

    helix_dictionary = {}

    i = -1
    helix_count = 0
    foundhelix = False
    for pos in list(bracket_notation):
        i += 1
        if pos == '(':
            lastRoundedPairedPos.append(i)
        elif pos == '{':
            lastCurlyPairedPos.append(i)
        elif pos == '[':
            lastSquarredPairedPos.append(i)

        # If this is found it means a new helix has occured
        elif pos == ')':
            if foundhelix is False:
                foundhelix = True
                helix_count += 1
                helix_dictionary[helix_count] = []
            temp_position = [lastRoundedPairedPos.pop(), i]
            temp_position.sort()
            helix_dictionary[helix_count].append(temp_position)
            # If the structure list is empty it means a helix has ened
            if len(lastRoundedPairedPos) == 0:
                foundhelix = False
        elif pos == '}':
            basePairs.append([lastCurlyPairedPos.pop(), i])
        elif pos == ']':
            basePairs.append([lastSquarredPairedPos.pop(), i])
        elif pos == '.':
            foundhelix = False

    # Now to find the cheaters
    additional_helices = []
    new_helix = []
    for helix in helix_dictionary:
        # Were going to see if these base indicies are contiguous
        first_round = True
        old_base0 = None
        old_base1 = None
        searching = False
        helix_copy = copy(helix_dictionary[helix])
        print helix_copy
        for base0, base1 in helix_copy:
            if searching is False:
                # Base0 is counting down base1 is counting up
                if first_round:
                    first_round = False
                    old_base0 = base0
                    old_base1 = base1
                    continue
                if (base0 != (old_base0 - 1)) or (base1 != (old_base1 + 1)):
                    # we found a cheater!
                    # We need to now search
                    searching = True
                    new_helix = []
                    # delete the one that was there
                    helix_dictionary[helix].remove([base0, base1])
                    new_helix.append([base0, base1])

            else:
                print base0, base1
                if (base0 != (old_base0 - 1)) or (base1 != (old_base1 + 1)):
                    # There's more than one cheater in this mix
                    additional_helices.append(new_helix)
                    # Do it again
                    new_helix = []
                    helix_dictionary[helix].remove([base0, base1])
                    new_helix.append([base0, base1])
                else:
                    helix_dictionary[helix].remove([base0, base1])
                    new_helix.append([base0, base1])
            # update our couters
            old_base0, old_base1 = base0, base1
        if len(new_helix) != 0:
            additional_helices.append(new_helix)
            new_helix = []
    # now to add the cheaters back as their own helixes
    for count, new_helix in enumerate(additional_helices):
        helix_dictionary[helix+1+count] = new_helix

    return helix_dictionary
