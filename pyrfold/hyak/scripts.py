"""
This module should contain all of the scripts which are needed
to build a set of python scripts that will be used to process data on kinefold
in a highthoughput fashion

Things it should do:
Delete node files
Delete other unneeded files
Convert .rnm files into timecourse files

"""
import os

def delete_nodes(nameofjobfile, nodedirectory, myscriptparmsdirectory):
    """
    :param nameofjobfile: name the job should have job will be added
    :type nameofjobfile: str
    :param nodedirectory: node directory path
    :type nodedirectory: str
    :param myscriptparmsdirectory: path to the myscript-parms directory for
        hyak submissions
    :type myscriptparmsdirectory: str:
    """
    if '.py' not in nameofjobfile:
        nameofjobfile += '.py'
    if 'job' not in nameofjobfile:
        nameofjobfile = 'job-' + nameofjobfile
    pathtofilename = os.path.join(myscriptparmsdirectory, nameofjobfile)
    with open(pathtofilename, 'wb') as pythonfile:
        modulesneeded = ['shutil']
        #write the needed modules
        for module in modulesneeded:
            pythonfile.write('import ' + module + '\n')
        #write the bulk of the work
        pythonfile.write('shutil.rmtree(' + nodedirectory + ')' + '\n')



