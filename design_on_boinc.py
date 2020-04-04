#!/usr/bin/env python
import os, sys
import string
import numpy as np
# from numpy import random
import re
# import pyrosetta
# from pyrosetta import *
# from pyrosetta.rosetta import *
import pandas as pd
import distutils.spawn
sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
import silent_tools
from silent_tools import eprint

from argparse import ArgumentParser

from collections import defaultdict
import time
import argparse
import itertools
import subprocess

import zipfile

from multiprocessing import Pool
import multiprocessing as mp


# Run these in bash to setup a working environment
#  source activate /software/conda/envs/tensorflow-testing (for numpy)
#  source /software/pyrosetta3.6/setup.sh  (for pyrosetta)
#
# This script generates blueprint, h-bond constraints, and Rosetta script xml files
# along with various utility scripts for running a local test, submitting jobs to
# R@h and getting results. This serves as an example and is not a production script.
#
# David K

# 
# Nate has edited this script so that it can take either a list of pdbs or a silent file, an xml file, flags file, number of jobs to be run, and pdbs per job
# It then writes this to a boinc submission script for the job specified by the pdbs and xml
# 
# Nate B

# 
# 
# Read in arguments and determine what kind of files we are working with
# 
# 

def cmd(command, wait=True, print_output=False):

    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    output = str(the_stuff[0]) + str(the_stuff[1])
    if ( print_output ):
        print(output)
    return output


if (len(sys.argv) == 1):
    eprint("")
    eprint('This script prepares your design jobs to run on Rosetta at Home')
    eprint("Arguments:")
    eprint(" -run_name          : optional, the name that you want to give to this whole job")
    eprint(" -xml               : the xml file to use, can be global or local path")
    eprint(" -flags             : the flags file to use, can be global or local path")
    eprint(" -pdbs_per_job      : optional, the number of structures to submit per boinc job")
    eprint(" -j                 : optional, how many cores you want this script to run on")
    eprint(" -in:file:silent    : the silent file of structures to design")
    eprint(" -extra_files       : space separated list of extra files to add to run")
    eprint(" -per_pdb_files     : score file format (tag last column) of extra files for each pdb. Use =-=> to rename file.")
    eprint(" -add_pdb_ids       : Add pdb ID to REMARK and script_vars")
    eprint(" -queue             : The queue setting for boinc. How many times should this job run?")
    eprint("                             or")
    eprint(" list of pdbs       : this script will also read in a list of pdbs if not given a silent file")
    sys.exit(1)

parser = ArgumentParser()
parser.add_argument("-run_name", type=str, default="design_on_boinc")
parser.add_argument("-xml", type=str, default="")
parser.add_argument("-flags", type=str, default="")
parser.add_argument("-pdbs_per_job", type=int, default=1)
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("-j", type=int, default=1)
parser.add_argument("-extra_files", type=str, nargs="*", default=[])
parser.add_argument("-per_pdb_files", type=str, default="")
parser.add_argument("-add_pdb_ids", action="store_true")
parser.add_argument("-priority", type=int, default=0)
parser.add_argument("-queue", type=int, default=1)
parser.add_argument('pdbs', type=str, help="input pdb", nargs="*")

args = parser.parse_args(sys.argv[1:])

run_filename = args.run_name
pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")
xml_filename = args.xml
flags_filename = args.flags
n_cores = args.j
pdbs_per_file = args.pdbs_per_job
per_pdb_files = args.per_pdb_files
add_pdb_ids = args.add_pdb_ids

if (xml_filename == '' or flags_filename == ''):
    sys.exit("This script needs both an xml file and a flags file")

using_silent = False
if (len(pdbs) == 0):
    using_silent = True
    if (silent == ''):
        sys.exit("This script needs either a list of pdb files or a silent file")

for file in args.extra_files:
    if ( not os.path.exists(file) ):
        sys.exit("%s doesn't exist"%file)
extra_files = args.extra_files

per_pdb_files_dict = defaultdict(list)
if ( per_pdb_files != "" ):
    with open(per_pdb_files) as f:
        for line in f:
            line = line.strip()
            if (len(line) == 0):
                continue
            sp = line.split()
            tag = sp[-1]

            # don't check for existence here as there might be a lot...
            for file in sp[:-1]:
                per_pdb_files_dict[tag].append(file)

if ( not os.path.exists("jobs") ):
    os.mkdir("jobs")
cmd("move_to_scratch jobs")

alpha = list(string.ascii_lowercase)
number = list(string.digits)

for num in number:
    for alph in alpha:
        fol = num + alph
        if ( not os.path.exists("jobs/%s"%fol) ):
            os.mkdir("jobs/%s"%fol)

# 
# 
# 
# Functions
# 
# 
# 


# This takes a list of pdbs and writes them to all_structs.silent
def deal_with_pdbs( list_o_pdbs ):
    pdbs_string = ' '.join( list_o_pdbs )
    cmd( '/home/nrbennet/software/silent_tools/silentfrompdbsparallel ' + pdbs_string + ' > all_structs.silent' )

def ensure_silent_is_binary( filename ):
    # could take this out and just read line in the file until I got to a new struct (would be faster)
    # Other than that this function works as intended
    # tags = cmd( '/home/nrbennet/software/silent_tools/silentls ' + filename, True )
    silent_index = silent_tools.get_silent_index( filename )
    tags = silent_index['tags']
    sample = tags[0]
    # if ( sample not in silent_index['index'] ):
    #         eprint("silentslice: Unable to find tag: %s"%sample)
    #         sys.exit("Something has gone seriously wrong")
    with open( filename ) as total_file:
        structure = silent_tools.get_silent_structure_file_open( total_file, silent_index, sample )
        for line in structure:
            if line[0] in 'ELH':
                return filename

    # cmd( 'touch binary.silent', False)
    # Do I need the extra res flag? -- no
    print("Converting non-binary silent to binary silent")

    command = 'score_jd2 -in:file:silent %s -out:file:silent binary.silent -out:file:silent_struct_type binary 1>&2'%filename

    cmd( command, True )
    return 'binary.silent'


class MyZip:

    def __init__(self, zip_name):
        self.zip_name = zip_name
        self.real_files = {}
        self.virtual_files = {}

    def add_real_file(self, fname, store_name=None):
        if ( store_name is None ):
            store_name = os.path.basename(fname)

        if ( store_name in self.real_files ):
            if (fname != self.real_files[store_name]):
                print("Error! zip name collision: %s from:"%(store_name))
                print("    "+real_file_info[store_name])
                print("    "+file)
            else:
                # here it's a true duplicate and we don't care
                return

        self.real_files[store_name] = fname

    def add_virtual_file(self, fname, contents):
        if ( fname in self.virtual_files ):
            if ( contents != self.virtual_files[fname]):
                print("Error! Name collision on virtual files: %s"%fname)

        self.virtual_files[fname] = contents



    def write(self):
        # print("MyZip: " + self.zip_name)
        with zipfile.ZipFile(self.zip_name, "w", zipfile.ZIP_DEFLATED) as z:
            for store_name in self.real_files:
                # print("    " + os.path.basename(real_file))
                real_file = self.real_files[store_name]
                z.write(real_file, store_name)
            for virtual_file in self.virtual_files:
                # print("    " + virtual_file)
                contents = self.virtual_files[virtual_file]
                z.writestr(virtual_file, contents)




input_scores_re = re.compile(r"\nREMARK[^\n]+_input_score")
id_re = re.compile(r"\nREMARK[^\n]+ID[^\n]+")

def make_run( xml_filename, runname, i, pdbs_per_file, silent_index, sf_open ):

    zip_name = "jobs/" + runname[:2] + "/" + runname + ".zip"
    myzip = MyZip(zip_name)

    myzip.add_real_file(xml_filename)
    for file in extra_files:
        myzip.add_real_file(file)

    extra_flags = ""
    id_vars = []
    extra_job_files = ""
    extra_job_files += ", " + zip_name

    new_silent_file = ""
    new_silent_file += silent_tools.silent_header(silent_index)

    structures = silent_tools.get_silent_structures_true_slice( sf_open, 
                silent_index, i, i+pdbs_per_file, True )
    tags = silent_index['tags'][i:i+pdbs_per_file]

    for i, pair in enumerate(zip(structures, tags)):
        s, t = pair

        first_n = s.find("\n")

        # this is a hack to allow variable length backbones on boinc
        #   specifically it tricks
        #   /projects/boinc/workspace/boinc/sched/validate_util.cpp
        #   on line 692 into skipping the coord_check
        if ( input_scores_re.search( s ) is None ):
            s = s[:first_n+1] + "REMARK _input_score 0\n" + s[first_n+1:]

        if ( add_pdb_ids ):
            s = id_re.sub("\n", s)
            s = s[:first_n+1] + "REMARK ID %s\n"%(t) + s[first_n+1:]
            id_vars.append("id%03i=%s"%(i, t))

        for item in per_pdb_files_dict[t]:
            if ( "=-=>" in item ):
                sp = item.split("=-=>")
                myzip.add_real_file(sp[0], sp[1])
            else:
                myzip.add_real_file(item)

        new_silent_file += s

    myzip.add_virtual_file('%s.silent'%runname, new_silent_file)


    if ( len(id_vars) > 0 ):
        temp = "jobs/%s/%s.flags"%(runname[:2], runname)
        extra_job_files += ", " + temp
        with open(temp, "w") as f:
            f.write("-script_vars " + " ".join(id_vars) + "\n")

        extra_flags += " @%s.flags"%runname


    return extra_flags, myzip, tags, extra_job_files



cores_available = os.getenv('SLURM_CPUS_ON_NODE') #cmd( 'echo $SLURM_CPUS_ON_NODE', True )
if cores_available is None:
    cores_available = mp.cpu_count()
else:
    cores_available = int( cores_available )

if ( n_cores > cores_available ):
    print( 'There are not %i cores available'%( n_cores ) )
    n_cores = cores_available
print( 'Using %i cores'%n_cores )

# this is how many commands have been run
init_count = mp.Value('i', 0)

def init_thread(counter_in):
    global counter
    counter = counter_in

def thread_write_zip(myzip):
    # global counter
    # counter.value += 1
    myzip.write()


# do this early to avoid sharing the silent_index
pool = Pool( n_cores, initializer=init_thread, initargs=(init_count,) )



# 
# 
# 
# Here is where it all gets put together and called
# 
# 
# 

assert( os.path.exists( flags_filename ) )
assert( os.path.exists( xml_filename ) )


total_flags_filename = os.path.realpath( flags_filename )

local_xml = os.path.basename(xml_filename) #cmd( 'basename %s'%xml_filename, True )
local_flag = os.path.basename(flags_filename) #cmd( 'basename %s'%flags_filename, True )


with open(xml_filename) as f:
    full_xml = f.read()
if ( "PoseComment" in full_xml and not args.add_pdb_ids):
    sys.exit("You probably forgot to specify -add_pdb_ids to commandline")

if ( "/home" in full_xml ):
    sys.exit("You specified a full path in your xml which is probably wrong")

if ( "~/" in full_xml ):
    sys.exit("\"~/\" appeared in your xml which is probably wrong")

# going to assume no spaces in paths so we can avoid flagging formulas
if ( not re.search(r'="[^" ]*/[^" ]+/', full_xml) is None ):
    sys.exit("Found a path in your xml which is probably wrong")
        

# total_xml_filename = os.path.realpath( xml_filename )
# scriptpath = os.path.dirname(os.path.realpath(__file__))

all_silent_filename = ''
if not using_silent:
    deal_with_pdbs( pdbs )
    all_silent_filename = all_structs.silent
else:
    all_silent_filename = silent
    # all_silent_filename = ensure_silent_is_binary( silent ) # my silent index takes 3 minutes to open...

print("Loading silent index")
silent_index = silent_tools.get_silent_index( all_silent_filename )

runnames = set()
tags = silent_index['tags'] #cmd( '/home/nrbennet/software/silent_tools/silentls ' + all_silent_filename, True )
num_tags = len(tags) #int(cmd( '/home/nrbennet/software/silent_tools/silentls %s | wc -l'%all_silent_filename, True ))


i = 0

run_extra_flags = {}
run_extra_files = {}

# doing it this way to avoid swear words
def random_name():
    # there has to be a better way...
    return ( np.random.choice(number, 1).item() + "".join(np.random.choice(alpha, 2))
            + np.random.choice(number, 1).item() + "".join(np.random.choice(alpha, 2))
            + np.random.choice(number, 1).item() + np.random.choice(alpha, 1).item()
            )

# This loop splits up tags in the order that they appear in the silent file
# it shamelessly steals all of that logic from bcov's silentsplitshuf
scheduled = 0
with open(all_silent_filename) as all_sf_open:
    with open(run_filename + ".info", "w") as f:
        while i < num_tags:

            # Don't let the main loop get too far ahead of the threads to keep memory low
            # if ( scheduled - init_count.value > n_cores * 20 ):
            #     time.sleep(0.1)
            #     continue

            if i + pdbs_per_file > num_tags:
                pdbs_per_file = num_tags - i
            runname = random_name()

            while ( runname in runnames ):
                runname = random_name()

            runname += "_" + run_filename

            runnames.add(runname)

            these_extra_flags, myzip, tags, these_extra_files = make_run( xml_filename, runname, i, pdbs_per_file, silent_index, all_sf_open )

            run_extra_flags[runname] = these_extra_flags
            run_extra_files[runname] = these_extra_files

            thread_write_zip(myzip)
            # pool.apply_async( thread_write_zip, args = ( myzip,) )
            scheduled += 1

            i += pdbs_per_file

            f.write(" ".join(tags) + " %s\n"%runname )


# This closes the pool and allows all the processes to complete before proceeding in
# the program
pool.close()
pool.join()


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


for ichunk, chunk in enumerate(chunks(list(runnames), 30000)):
    # Here is where I actually write the boinc file
    with open( run_filename + '_%i.boinc'%ichunk , 'w' ) as f:
        for i in chunk:
            f.write("application = rosetta\npriority = %i\n\n"%args.priority)
            # f.write("name = %s\n"%run_filename)
            f.write("name = %s_SAVE_ALL_OUT\n"%i)
            f.write("description = %s\n"%i)
            f.write("inputfiles = %s%s\n"%(total_flags_filename, run_extra_files[i]))
            f.write("arguments = -run:protocol jd2_scripting -parser:protocol " + local_xml + 
                # " -out:prefix " + i + "_" +
                " -database minirosetta_database @" + 
                local_flag + " -in:file:silent " + i + ".silent -in:file:silent_struct_type binary -silent_gz -mute all" +
                " -out:file:silent_struct_type binary -out:file:silent default.out -in:file:boinc_wu_zip " +
                i + ".zip" + 
                run_extra_flags[i] + 
                "\n")
            
            f.write("resultfiles = default.out.gz\n")
            # May want to change this later but for now I'm leaving it hard-coded as 1
            f.write("queue = %i\n"%args.queue)


num_chunks = ichunk+1

# 
# 
# 
# Writing convenience scripts
# 
# 
# 


# create test shell script
f = open( run_filename + '.test', 'w')
f.write('/projects/boinc/bin/run_test_rah.pl ' + run_filename + '_0.boinc')
f.close()

# create boinc submit shell script
f = open( run_filename + '.submit', 'w')
for ichunk in range(num_chunks):
    f.write("/projects/boinc/bin/boinc_submit " + run_filename + "_%i.boinc\n"%(ichunk))
f.write("\n# Useful scripts (run with - for usage)\n")
f.write("#   /projects/boinc/bin/boinc_q <batchid>\n")
f.write("#   /projects/boinc/bin/boinc_resize -size <size> <batchid>\n")
f.write("#   /projects/boinc/bin/boinc_rm <batchid>\n")
f.write("# Status page\n")
f.write("#   https://boinc.bakerlab.org/queue (username: boinc password: results)\n")
f.close()

# create results retrieval shell script
# probably will need to modify this -NRB

# f = open( run_filename + '.results', 'w')
# f.write("cp /projects/boinc-results/"+job_name[:8]+"/"+job_name+"*.all.out.bz2 ./\n")
# f.write("bunzip2 "+job_name+"*.all.out.bz2\n")
# f.write("\n# Extract silent file for pdbs\n/software/rosetta/latest/bin/extract_pdbs -in:file:silent_struct_type binary -silent_read_through_errors -in:file:silent "+job_name+"*.all.out\n")
# f.write("\n# Here is an example of running various design selection filters. This will also extract the pdbs.\n");
# f.write("/home/dekim/src/design_on_boinc/Rosetta_dev/main/source/bin/rosetta_scripts.hdf5.linuxgccrelease -in:file:silent "+job_name+"*.all.out -parser:protocol 
#/home/dekim/src/design_on_boinc/inputfiles/filter.xml -beta @/home/dekim/src/design_on_boinc/inputfiles/filter_flags -out:file:scorefile "+job_name+".sc\n")
# f.close()



