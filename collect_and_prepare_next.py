#!/usr/bin/env python

import sys
import os
import getpass

import argparse
from collections import defaultdict

import re


import distutils.spawn
sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
import silent_tools
from silent_tools import cmd
import random


parser = argparse.ArgumentParser()
parser.add_argument("-input_silents", type=str, nargs="+")
# parser.add_argument("-batch_id_map", type=bool, nargs="+")
parser.add_argument("-info_files", type=str, nargs="+")
parser.add_argument("-collected_jobs", type=str, nargs="*", default=[])
parser.add_argument("-collected_silents", type=str, nargs="*", default=[])
parser.add_argument("-resubmit_pending_tags", action="store_true")
parser.add_argument("-collect_partially_complete", action="store_true")
parser.add_argument("-assume_missing_means_done", action="store_true")
parser.add_argument("-force_recollect_patterns", type=str, nargs="*", default=[])
parser.add_argument("-recollect_verbose_frac", type=float, default=0.001)
parser.add_argument("-user", type=str, default="")
parser.add_argument("-j", type=int, default=1)

args = parser.parse_args(sys.argv[1:])


user = args.user
if ( user == "" ):
    user = getpass.getuser()



cmd("touch scratch_me")
cmd("move_to_scratch scratch_me")
scratch_path = os.path.dirname(cmd("readlink -f scratch_me")) + "/"

print("Using scratch path: %s"%scratch_path)



all_input_tags = set()

if ( hasattr( os, "all_input_tags" ) ):
    all_input_tags = os.all_input_tags

else:
    print("Reading input tags")
    for input_silent in args.input_silents:
        print("    " + input_silent)
        silent_index = silent_tools.get_silent_index( input_silent )

        for tag in silent_index['tags']:
            all_input_tags.add(tag)
    os.all_input_tags = all_input_tags

# print("Reading batch ids")
# id_to_job = {}
# with open(args.batch_id_map) as f:
#     for line in f:
#         line = line.strip()
#         if (len(line) == 0):
#             continue
#         sp = line.split()

#         for idd in sp[:-1]:
#             id_to_job[idd] = sp[-1]

print("Reading info files")
tag_to_runnames = defaultdict(list)
collectable_runnames = set()

if ( hasattr( os, "tag_to_runnames" ) ):
    tag_to_runnames = os.tag_to_runnames
else:
    for info_file in args.info_files:
        print("    " + info_file)

        job = os.path.basename(info_file)
        assert( job.endswith(".info") )
        job = job[:-len(".info")]

        with open(info_file) as f:
            for line in f:
                line = line.strip()
                if ( len(line) == 0):
                    continue
                sp = line.split()

                runname = sp[-1]
                collectable_runnames.add(runname)

                for tag in sp[:-1]:
                    tag_to_runnames[tag].append(runname)
    os.tag_to_runnames = tag_to_runnames

recollect_patterns = []
for pattern in args.force_recollect_patterns:
    recollect_patterns.append(re.compile(pattern))

print("Reading collected jobs")
collected_runnames = set()
for collected_job_file in args.collected_jobs:
    print(collected_job_file)
    with open(collected_job_file) as f:
        for line in f:
            line = line.strip()
            if (len(line) == 0):
                continue
            for pattern in recollect_patterns:
                if ( not pattern.search(line) is None ):
                    if ( random.random() < args.recollect_verbose_frac ):
                        print("Recollecting: %s"%line)
                    continue
            collected_runnames.add(line)

def parse_silent_file(fname, tags_set):
    print("    " + fname)
    silent_index = silent_tools.get_silent_index( fname )

    for tag in silent_index['orig_tags']:
        tag = tag[:-len("00004_0000001_0")]
        tags_set.add(tag)

print("Reading collected tags")
collected_tags = set()
for input_silent in args.collected_silents:
    parse_silent_file(input_silent, collected_tags)




print("Getting boinc status")
if ( hasattr(os, "boinc_q_raw") ):
    boinc_q_raw = os.boinc_q_raw
else:
    boinc_q_raw = cmd("/home/bcov/util/boinc_q_bcov.py -long_names")
    os.boinc_q_raw = boinc_q_raw

boinc_q_raw = boinc_q_raw.strip().split("\n")
headers = boinc_q_raw[0].split()
to_col = {}
for i, head in enumerate(headers):
    to_col[head] = i

runname_to_status = {}
runname_to_fullname = {}
runname_to_batch = {}
for line in boinc_q_raw[1:]:
    line = line.strip()
    if (len(line) == 0):
        continue
    sp = line.split()
    status = sp[to_col["S"]]
    runname = sp[to_col["NAME"]]

    runname = runname.replace("_SAVE_ALL_OUT", "")

    runname_to_fullname[runname] = sp[to_col["NAME"]]
    runname_to_batch[runname] = sp[to_col["ID"]].split(".")[0]

    if ( status == "X" ):
        status = "C"

    if ( int(sp[to_col["DECOYS"]]) == 0 ):
        status = "N"
    else:
        if ( args.collect_partially_complete ):
            status = "C"

    runname_to_status[runname] = status


# now we work on derived products

to_collect = []
for runname in collectable_runnames:
    if ( runname in collected_runnames ):
        continue
    if ( runname not in runname_to_status ):
        print("Missing job?: %s"%runname)
        if ( args.assume_missing_means_done ):
            to_collect.append(runname)
            runname_to_fullname[runname] = runname + "_SAVE_ALL_OUT"
            runname_to_batch[runname] = "*"
        continue
    if ( runname_to_status[runname] == "C" ):
        to_collect.append(runname)

# to_collect = to_collect[:10]

all_paths = []
for runname in to_collect:
    fullname = runname_to_fullname[runname]
    prefix = fullname[:8]
    batch = runname_to_batch[runname]

    has_all = ".all" if fullname != runname else ""

    path = "/projects/boinc-results/%s/%s_%s_0%s.out.bz2"%(
                prefix, fullname, batch, has_all)
    all_paths.append(path)

print("Collecting %i results"%len(all_paths))

with open("tmp", "w") as f:
    for path in all_paths:
        f.write("%s\n"%(path))

print(cmd("cat tmp | parallel -j %i --xargs 'eval bzcat' > %stmp.silent"%(args.j, scratch_path)))
os.remove("tmp")

print("Removing corrupt models")
print(cmd("silentdropcorruptmodels %stmp.silent > %stmp2.silent"%(scratch_path, scratch_path)))
os.remove("%stmp.silent"%scratch_path)
os.remove("%stmp.silent.idx"%scratch_path)

print("Renaming duplicate names")
print(cmd("silentls %stmp2.silent | silentrename %stmp2.silent > %scollected.silent"%(scratch_path, scratch_path, scratch_path)))
os.remove("%stmp2.silent"%scratch_path)
os.remove("%stmp2.silent.idx"%scratch_path)

cmd("ln -s %scollected.silent ."%(scratch_path))


print("Reading new tags")
newly_collected_tags = set()
parse_silent_file("collected.silent", newly_collected_tags)

print("Found %i new tags"%(len(newly_collected_tags - collected_tags)))

print("Determining tags that need to be resubmitted")
missing_tags = all_input_tags - newly_collected_tags - collected_tags

if ( args.resubmit_pending_tags ):
    tags_to_submit = missing_tags
else:
    tags_to_submit = set()
    for tag in missing_tags:
        runnames = tag_to_runnames[tag]
        is_pending = False
        for runname in runnames:
            if ( runname not in runname_to_status ):
                # print("Missing job?: %s")
                continue
            if ( runname_to_status[runname] in "RI" ):
                is_pending = True

        if ( not is_pending ):
            tags_to_submit.add(tag)


with open("to_resubmit.list", "w") as f:
    for tag in tags_to_submit:
        f.write("%s\n"%tag)

with open("collected.list", "w") as f:
    for runname in to_collect:
        f.write("%s\n"%runname)

print("Need to submit %i more tags"%len(tags_to_submit))




