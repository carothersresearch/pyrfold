"""module loads picked data to be used
RNABINDINGDICT is [sizeofmer][sequence][list, of, sequences, which, bind]

As of now all of the data that is presented in this init file was calcualted
using ViennaRNA RNAduplex

RNAPROBBINDING_345 : estimates a relative chance of binding by consiering
    the thermodynamics of each mers binding and putting it through a general
    connonical function

RNAINTERACTION : lists all of the interactions that are possible between given
    mers. This considers all wobble pairigs A:U G:C G:U

All of the calculations made here do not assume anything about single base
mismatches nor thermodynamics
"""
import cPickle as pickle
import os

with open(os.path.join(os.path.dirname(__file__),
                      '3&4&5mer_canonical_probability.p'), 'rb') as temppickle:
    RNAPROBBINDING_345 = pickle.load(temppickle)

with open(os.path.join(os.path.dirname(__file__),
                        '3&4mer_canonical_probability.p'), 'rb') as temppickle:
    RNAPROBBINDING_34 = pickle.load(temppickle)

with open(os.path.join(os.path.dirname(__file__),
                        '4&5mer_canonical_probability.p'), 'rb') as temppickle:
    RNAPROBBINDING_45 = pickle.load(temppickle)

with open(os.path.join(os.path.dirname(__file__),
                        '4&5mer_canonical_probability.p'), 'rb') as temppickle:
    RNAPROBBINDING_45 = pickle.load(temppickle)

with open(os.path.join(os.path.dirname(__file__),
                   '3&4&5mer_canonical_thermo_interac.p'), 'rb') as temppickle:
    RNATHERMO_345 = pickle.load(temppickle)

with open(os.path.join(os.path.dirname(__file__),
                     '3&4&5mer_canonical_interactions.p'), 'rb') as temppickle:
    RNAINTERACTION = pickle.load(temppickle)
