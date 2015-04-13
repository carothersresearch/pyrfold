"""
This package contains a wrapper for the kinefold binary
"""
import os
import subprocess
from tempfile import mkdtemp
from os.path import isdir
from shutil import rmtree
from ..hyak import process
import random

devnull = open(os.devnull, 'w')


#First the location of the kinefold biany needs to be found
path_to_kinefold = os.path.join(os.path.dirname((os.path.realpath(__file__))),
                                'applications',
                                'kinefold_long_static')

class Kinefold():
    """
    Object for calculating kinefold folding for RNA
    :param sequence: str
    :type sequence: string of RNA to be processed
    :param polrate: the polrate
    :type polrate: float
    :param dwell_time: dwell time after the polymerase has elongated
    :type dwell_time: float
    """
    def __init__(self, FoldingSubData):
        self._tempdir = ''
        self.FoldingSubData = FoldingSubData
        self.name = FoldingSubData.name
        self.timecoursedata = None
        self.compressedtimecoursedata = None

    def run_simulations(self):
        # confirm that there is a temp directory to write the files to
        self._check_tempdir()
        # Make the .dat file
        write_dat_files(self._tempdir, [self.name],
                        {self.name:self.FoldingSubData})
        write_req_files(self._tempdir, self._tempdir, self._tempdir,
                        [self.name], {self.name:self.FoldingSubData},
                        wrapper_run=True)
        # list all of the req files
        list_of_req_files = [os.path.join(self._tempdir, fi) for
                             fi in os.listdir(self._tempdir)
                             if (os.path.isfile(os.path.join(self._tempdir, fi))
                                 and '.req' in fi)]
        for reqfile in list_of_req_files:
            subprocess.call([path_to_kinefold, reqfile, '-noprint'],
                            stdout=devnull, stderr=devnull)
        # Process the data
        self.timecoursedata, self.compressedtimecoursedata = \
            process.timecourse(self._tempdir, 'dummy', singledirectory=True,
                               return_dictionary=True)
        # delete the directory
        self._close()

    def _close(self):
        '''Close the temporary dir (keeps /tmp clean).'''
        rmtree(self._tempdir)

    def _check_tempdir(self):
        '''If temp dir has been removed, create a new one.'''
        if not isdir(self._tempdir):
            self._tempdir = mkdtemp()


def dna_to_rna(seq):
    """(str) -> changed string
    simple function to replace all T with U
    """
    seq = seq.upper()
    seq = seq.replace("T","U")
    return seq


def write_dat_files(datdirectory, listofdevices, devicenametosubobj,
                    return_path=False):
    """This will fill in the tempnode/dat directory with the neccessary
    sequence information for processing"""
    for device in listofdevices:
        with open(os.path.join(datdirectory, (device + '.dat')), 'wb') as f:
            f.write("<" + device + '\n')
            # truncating device from tuple (seq,winStart,winStop)
            seq = devicenametosubobj[device].sequence
            start = devicenametosubobj[device].windowstart
            stop = devicenametosubobj[device].windowstop
            # accounting for shift in frame str[n:c] doesn't actually read
            # through to c it will stop at c-1
            seq = seq[start - 1: stop]
            seq = dna_to_rna(seq)
            f.write(seq + '\n')
    if return_path:
        return os.path.join(datdirectory, (device + '.dat'))


def write_req_files(parmdirectory, outputpath, datdirectory, listofdevices,
                    devicenametosubobj, wrapper_run=False):
    """(str, int, tuple) -> write .req files for device
    This will write all of the neccessary req files for kinefold to process
    """
    random.seed()  # uses system time to initialize random generator
    # EXPERIMENTAL VARIABLES
    for devicename in listofdevices:
        numberofsimulations = devicenametosubobj[devicename].numberofsimulations
        polrate, requestedtime = devicenametosubobj[devicename].kine_folding_data()
        psudoknots = devicenametosubobj[devicename].pseudoknots
        entanglements = devicenametosubobj[devicename].entanglements
        forcedhelixes = devicenametosubobj[devicename].forcedhelixes
        exptype = devicenametosubobj[devicename].experimenttype
        helix_min_free_eng = devicenametosubobj[devicename].helixminfreeenergy
        for i in range(numberofsimulations):
            # create name for req files and outputs
            reqname = str(i + 1).zfill(3) + devicename
            # open unique req file
            with open(os.path.join(parmdirectory, 'job.' + reqname
                      + '.req'), 'wb') as f:
                # Write random number seed
                f.write(str(int(round(random.random() * 10000))).zfill(4))
                # Write all output directories
                fileext = ['.p', '.e', '.rnm', '.rnms', '.rnml', '.rnm2']
                for ext in fileext:
                    f.write('\n')
                    if wrapper_run:
                        f.write(os.path.join(outputpath, reqname + ext))
                    else:
                        f.write(os.path.join(outputpath, devicename,
                                             reqname + ext))
                # Write .dat directory
                f.write('\n' + os.path.join(datdirectory, devicename + '.dat'))
                # Write 0 RNA 1 for DNA
                f.write('\n' + str(0))
                # helix minimum free energy in kcal/mol: 6.3460741=10kT
                f.write('\n' + str(helix_min_free_eng))
                # just something needed
                f.write('\n' + str(10000000))
                # requested folding time
                f.write('\n' + str(requestedtime))
                # psudoknots 1 = yes, 0 = no
                f.write('\n' + str(int(psudoknots)))
                # entanglements 1 = yes, 0 = no
                f.write('\n' + str(int(entanglements)))
                # simulation type: 1=renaturation; 2 20 =cotrans. @ 20msec/nt
                if exptype == 2:
                    f.write('\n' + '2 ' + str(polrate))
                else:
                    f.write('\n' + '1')
                for forc in forcedhelixes:
                    f.write('\n' + 'F ' + str(forc[0]) + ' ' + str(forc[1]) +
                                                         ' ' + str(forc[2]))
                #filename and filename.zip
                f.write('\n' + devicename + '\n' + devicename + '.zip' + '\n')


def run_simulations(path_to_req_files, name_of_device):
        # confirm that there is a temp directory to write the files to
        list_of_req_files = [os.path.join(path_to_req_files, fi) for
                             fi in os.listdir(path_to_req_files)
                             if (os.path.isfile(os.path.join(path_to_req_files, fi))
                                 and name_of_device in fi and '.req' in fi)]
        for reqfile in list_of_req_files:
            subprocess.call([path_to_kinefold, reqfile, '-noprint'],
                            stdout=devnull, stderr=devnull)
        # Delete the req files
        for reqfile in list_of_req_files:
            os.remove(reqfile)
