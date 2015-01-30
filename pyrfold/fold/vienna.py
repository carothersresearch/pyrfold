'''Vienna RNA module.'''
from subprocess import Popen, PIPE, STDOUT
from tempfile import mkdtemp
from os.path import isdir
from shutil import rmtree


class Vienna(object):
    '''Run Vienna RNA functions on a sequence.'''
    def __init__(self, seqs, dotbrackets=None, constraintstructures=None):
        '''
        :param seq: DNA or RNA sequences to evaluate.
        :type seq: list of pymbt.DNA, pymbt.RNA, or str
        :param dotbrackets: Dot bracket formatted structures of given seuqneces
        :type dotbrackets: list of str

        :returns: pymbt.analysis.Vienna instance.

        '''
        self._seqs = [str(seq) for seq in seqs]
        if dotbrackets:
            self._dotbrackets = dotbrackets
        else:
            self._dotbrackets = []
        if constraintstructures:
            self.constraintstructures = constraintstructures
        else:
            self.constraintstructures = []

    def free_energy_of_structure(self, temp=37, index=None):
        """Calculate the free energy of a given sequence and structure
        at a given temperature
        :param temp: Temperature at which to run calculations.
        :type temp: float
        :param index: Specific sequence index on which to calcualte mfe.
        :type index: int
        :returns: Minimum Free Energy (mfe) list [kcal/mol].
        :rtype: list of floats
        """
        if index is None:
            tempseqs = self._seqs
            tempdotbracets = self._dotbrackets
        else:
            tempseqs = [self._seqs[index]]
            tempdotbracets = [self._dotbrackets[index]]
        # self._check_tempdir()
        mfes = []
        for counter, seq in enumerate(tempseqs):
            process = Popen(['RNAeval', '-T', str(temp)], stdin=PIPE,
                            stdout=PIPE, stderr=STDOUT, cwd=self._tempdir)
            tempinput = seq + '\n' + tempdotbracets[counter]
            output = process.communicate(input=tempinput)[0]
            lines = output.splitlines()
            lines = lines[-1].split('(')[-1].split(')')[0].strip()
            mfes.append(float(lines))
        # self._close()
        return mfes

    def mfe(self, temp=50.0, returnstructure=False, index=None,
            constraintstructure=False):
        '''Calculate the minimum free energy.
        :param temp: Temperature at which to run calculations.
        :type temp: float
        :param index: Specific sequence index on which to calcualte mfe.
        :type index: int
        :returns: Minimum Free Energy (mfe) list [kcal/mol].
        :rtype: list of floats

        '''
        mfes = []
        dotbrackets = []
        ARGUMENT = ['RNAfold', '-T', str(temp), '--noPS']
        if constraintstructure:
            ARGUMENT.append('-C')
            if index is None:
                tempseqs = self._seqs
                tempconstaints = self.constraintstructures
            else:
                tempseqs = [self._seqs[index]]
                tempconstaints = [self.constraintstructures[index]]
            for seq, const in zip(tempseqs, tempconstaints):
                process = Popen(ARGUMENT, stdin=PIPE,
                                stdout=PIPE, stderr=STDOUT)
                stringtoinput = seq + '\n' + const
                output = process.communicate(input=stringtoinput)[0]
                lines = output.splitlines()
                mfe = lines[-1].split('(')[-1].split(')')[0].strip()
                dotbracket = lines[-1].split(' ')[0].strip()
                dotbrackets.append(dotbracket)
                mfes.append(float(mfe))
        else:
            if index is None:
                tempseqs = self._seqs
            else:
                tempseqs = [self._seqs[index]]
            # self._check_tempdir()
            for seq in tempseqs:
                process = Popen(ARGUMENT, stdin=PIPE,
                                stdout=PIPE, stderr=STDOUT)
                output = process.communicate(input=seq)[0]
                lines = output.splitlines()
                mfe = lines[-1].split('(')[-1].split(')')[0].strip()
                dotbracket = lines[-1].split(' ')[0].strip()
                dotbrackets.append(dotbracket)
                mfes.append(float(mfe))
            # self._close()
        if returnstructure:
            return zip(mfes, dotbrackets)
        return mfes

    def pairs(self, temp=50.0, index=None):
        '''Calculate per-pair probability of being unbound (secondary
        structure).
        :param temp: Temperature at which to run calculations.
        :type temp: float
        :param index: Specific sequence index on which to calcualte pairs.
        :type index: int
        :returns: Pair probability for every base in the sequences.
        :rtype: list of list of floats.

        '''
        if index is None:
            tempseqs = self._seqs
        else:
            tempseqs = [self._seqs[index]]
        self._check_tempdir()
        unbounds = []
        for seq in tempseqs:
            process = Popen(['RNAfold', '-p', '-T', str(temp)], stdin=PIPE,
                            stdout=PIPE, stderr=STDOUT, cwd=self._tempdir)
            process.communicate(input=seq)[0]
            with open('{}/dot.ps'.format(self._tempdir), 'r') as dot:
                text = dot.read()
                text = text[text.index('%data'):]
            split = text.split('\n')
            data = [x for x in split if x.endswith('ubox')]
            data = [x.rstrip(' ubox') for x in data]
            data = [x.split() for x in data]
            data = [(int(a), int(b), float(c)) for a, b, c in data]
            unbound = [1.0] * len(seq)
            for base1, base2, prob_sqr in data:
                probability = prob_sqr**2
                unbound[base1 - 1] -= probability
                unbound[base2 - 1] -= probability
            unbounds.append(unbound)
        self._close()
        return unbounds

    def hybridization_mfe(self, simtype='duplex', bindingsequence=False,
                 temp=37.0, indexs=None, noGU=False):
        """ Calculate the mfe of two RNA seqeunces complexing
        :param temp: Temperature at which to run calculations.
        :type temp: float
        :param indexs: list (len 2) of Specific sequence index on which to
            calcualte pairs.
        :type index: list
        :param simtype: Vienna simulation package to use. 'duplex' dose not allow
            self RNA interactions. 'cofold' allows self RNA interactions to
            occur
        :type simtype: str
        :param bindingsequence: If true output will be a touple contining
        (MFE, [dotbracket1, dotbracket2] [[start1,stop1], [[start2,stop2]])
            Where 1 and 2 are the first and second indexs processed
        :type bindingsequence: bool
        :returns: Binding energy between two given sequence
        :rtype: list of list of floats.

        TODO: Write bindingsequence for cofold
        """
        if indexs is None:
            tempseqs = self._seqs
        else:
            tempseqs = [self._seqs[indexs[0]], self._seqs[indexs[1]]]
        # self._check_tempdir()
        sub_command = []
        if simtype == 'cofold':
            sub_command.append('RNAcofold')
        elif simtype == 'duplex':
            sub_command.append('RNAduplex')
        else:
            #Write some error thing here
            pass
        sub_command.append('-T')
        sub_command.append(str(temp))
        if noGU:
            sub_command.append('--noGU')
        if simtype == 'cofold':
            process = Popen(sub_command, stdin=PIPE,
                            stdout=PIPE, stderr=STDOUT)
            stringtoinput = tempseqs[0] + '&' + tempseqs[1]
            output = process.communicate(input=stringtoinput)[0]
            lines = output.splitlines()
            mfe = float(lines[-1].split('(')[-1].split(')')[0].strip())
        elif simtype == 'duplex':
            process = Popen(sub_command, stdin=PIPE,
                            stdout=PIPE, stderr=STDOUT)
            stringtoinput = tempseqs[0] + '\n' + tempseqs[1]
            output = process.communicate(input=stringtoinput)[0]
            mfe = float(output.split('(')[-1].split(')')[0])
            brackets = [output.split('&')[0], output.split('&')[1].split(' ')[0]]
            bindlocations = [[int(i)
                         for i in output.split('  ')[1].strip().split(',')],
                        [int(i)
                         for i in output.split('  ')[3].strip().split(',')]]
        # self._close()
        if bindingsequence:
            return [(mfe, brackets, bindlocations)]
        return [[mfe]]

    def _close(self):
        '''Close the temporary dir (keeps /tmp clean).'''
        rmtree(self._tempdir)

    def _check_tempdir(self):
        '''If temp dir has been removed, create a new one.'''
        if not isdir(self._tempdir):
            self._tempdir = mkdtemp()
