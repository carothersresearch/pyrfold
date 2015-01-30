"""
This is a script which will be called many times (as many times as there are
unique devices that are submitted to hyak)

This script will parse the submission file for sequence, polrate, number of
simulations, sequence to fold, forced helicies (if any)
"""
from sys import argv
import cPickle as pickle
import os
import pyrfold.hyak.process as hyakp
import pyrfold.fold.kinefold as kinefold
import pyrfold.pyrfile as pyrfile
import shutil

# This script is expecting a pickle file to be sent to it
pickle_file_path = argv[1]
if 'final_compress' in pickle_file_path:
    tempfilepath = argv[2]
    outputfilepath = argv[3]
    processoutputpath = argv[4]

    timecoursepath = os.path.join(processoutputpath, 'timecourse')
    compresseddictpath = os.path.join(processoutputpath, 'compressedtime')

    NODEPATH = os.path.split(os.path.abspath(tempfilepath))[0]
    shutil.rmtree(NODEPATH)

    # Check to see if there are any more node fi
    exppath = os.path.split(os.path.abspath(NODEPATH))[0]
    nodelist = [ls for ls in os.listdir(exppath)
                if os.path.isdir(os.path.join(exppath, ls)) and
                'node' in ls]

    if len(nodelist) == 0:
        # If there are no ondefiles - it is safe to compress everything
        hyakp.compress_and_delete_directory(outputfilepath)
        hyakp.compress_and_delete_directory(timecoursepath)
        hyakp.compress_and_delete_directory(compresseddictpath)

else:
    with open(pickle_file_path, 'rb') as temppick:
        kinesubobj = pickle.load(temppick)

    # Make a dummy dict so everything is compliant with old code
    name = kinesubobj.name
    nametokinesubobj = {name:
                        kinesubobj}

    tempfilepath = argv[2]
    outputfilepath = argv[3]
    processoutputpath = argv[4]

    # Let's grab the sub summary
    SUBSUMMARYFILE = os.path.split(os.path.abspath(processoutputpath))[0]
    SUBSUMMARYFILE = os.path.join(SUBSUMMARYFILE, 'sub_summary.p')
    SUBSUMMARYDATA = pyrfile.load_pickled_sub_summary(SUBSUMMARYFILE)
    # Now that we have the kinsuboj file we can do work on it


    # STEP 1 Write the dat file
    dat_path = \
        kinefold.write_dat_files(tempfilepath, [name], nametokinesubobj,
                                 return_path=True)

    # STEP 2 make directory for all simulation data (in the output folder)
    rawdataoutput = os.path.join(outputfilepath, name)
    os.mkdir(rawdataoutput)

    # STEP 3 Make and run all require REQ files (in the output folder)
    kinefold.write_req_files(tempfilepath, outputfilepath, tempfilepath, [name],
                             nametokinesubobj, wrapper_run=False)
    kinefold.run_simulations(tempfilepath, name)

    # STEP 4 Delete all unneeded data (in the output folder and temp)
    hyakp.clean_up_output_data(rawdataoutput, singlefile=True)
    hyakp.delete_files_in_folder(tempfilepath, name)

    # STEP 5 Process the data to create (Transfer to proc_output)
    hyakp.timecourse(rawdataoutput, processoutputpath, singledirectory=True)

    # STEP 6 Compress output/device data
    hyakp.compress_and_delete_directory(rawdataoutput)


