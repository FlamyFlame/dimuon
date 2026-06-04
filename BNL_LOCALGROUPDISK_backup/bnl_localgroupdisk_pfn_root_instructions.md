# BNL LOCALGROUPDISK PFN/ROOT Notes From Luke

Source: Mattermost DM with Luke Jarrett Mozarsky on May 19. I accessed the visible Mattermost thread and downloaded the five attachments that were present in the May 19 exchange. The Mattermost day divider showed `May 19`; the app accessibility text rendered the messages as `Tuesday, May 19`.

Local attachment copies:

- [pfns_ROOT_quickstart_guide.pdf](mattermost_luke_may19_attachments/pfns_ROOT_quickstart_guide.pdf)
- [generateArgs.sh](mattermost_luke_may19_attachments/generateArgs.sh)
- [pfns.txt](mattermost_luke_may19_attachments/pfns.txt)
- [condor.sub](mattermost_luke_may19_attachments/condor.sub)
- [run.sh](mattermost_luke_may19_attachments/run.sh)

## May 19 Chat Summary

Yuhan asked Luke two things:

- Whether Luke had slides or sample analysis code for running analysis directly on files in the BNL local group disk area using pointers in dCache/PNFS.
- Whether Luke had CC and PC streams skimmed for any Run 3 PbPb year with muon trigger decisions saved.

Luke replied that he had an example macro/workflow for reading files using PFNs and would send it. For the CC/PC stream question, he said he did not think he had those skims with muon triggers, but would double check.

Later that day Luke sent the PFN/ROOT workflow. His message said the first step is to make sure the dataset is stored on disk on the RSE `BNL-OSG2_LOCALGROUPDISK`. Then, after setting up Rucio, write each file PFN into a text file with:

```bash
rucio list-file-replicas <scope:DID> --protocols root --pfns --rses BNL-OSG2_LOCALGROUPDISK > <file_name>.txt
```

Luke then said what to do next depends on whether Condor is being used:

- Without Condor: read one line at a time from the PFN text file and add each PFN to a `TChain` as with locally stored files.
- With Condor: likely submit one job per file, with additional per-job arguments as needed.

Luke's Condor example used `generateArgs.sh` to convert `pfns.txt` into `args.txt`, then `condor.sub` to submit jobs that execute `run.sh`. He noted two points for Condor:

- Make sure the grid proxy stays valid for the duration of all Condor jobs.
- Ensure the correct grid proxy location is passed in `transfer_input_files` in `condor.sub`.

## PDF Attachment Notes

Attachment: [pfns_ROOT_quickstart_guide.pdf](mattermost_luke_may19_attachments/pfns_ROOT_quickstart_guide.pdf)

The PDF title page says:

```text
Quickstart guide to using rucio PFNS with ROOT
March 31, 2026
Luke Mozarsky
Columbia HI Weekly Meeting
```

The PDF step slide says:

1. For work on BNL, ensure the dataset is stored on the `BNL-OSG2_LOCALGROUPDISK` storage element.
2. This is "not strictly necessary," but the slide says speeds may slow down significantly if datasets are stored on a non-BNL RSE.
3. Make a text file containing PFNs with:

```bash
rucio list-file-replicas <scope:DID> --protocols root --pfns --rses BNL-OSG2_LOCALGROUPDISK > <file_name>.txt
```

The slide also says it is optional to add this to `.bashrc` or `.zshrc`:

```bash
export SITE_NAME=BNL-ATLAS
```

The slide attributes that recommendation to Mason Housenga and says it may speed up the Rucio download process, but Luke had not personally tested it.

The PDF states that from this point, paths in `<file_name>.txt` can be treated as if they were local paths, for example if they contain nTuples, they can be added to a `TChain` and manipulated as usual.

The PDF also states that ROOT needs Rucio set up with a valid grid proxy for the full runtime. It gives the example that if code will run longer than 12 hours, the default `voms-proxy-init` lifetime may need manual adjustment.

## Non-Condor Usage

This section follows Luke's chat message and the PDF attachment.

1. Set up Rucio and a valid ATLAS grid proxy.
2. Confirm the dataset has replicas on `BNL-OSG2_LOCALGROUPDISK`.
3. Produce a PFN text file with `rucio list-file-replicas`.
4. In ROOT code, read the PFN text file line-by-line and add each PFN to the `TChain`.

No standalone non-Condor macro was attached in the May 19 thread. The `TChain` usage above is a summary of Luke's message and PDF wording, not copied from an attached macro.

## Condor Usage From Attachments

The Condor workflow is represented by four downloaded attachments:

- [pfns.txt](mattermost_luke_may19_attachments/pfns.txt)
- [generateArgs.sh](mattermost_luke_may19_attachments/generateArgs.sh)
- [condor.sub](mattermost_luke_may19_attachments/condor.sub)
- [run.sh](mattermost_luke_may19_attachments/run.sh)

The intended flow from Luke's chat is:

```bash
sh generateArgs.sh
condor_submit condor.sub
```

The attachment `generateArgs.sh` creates `args.txt` from `pfns.txt`. The attachment `condor.sub` reads rows from `args.txt` and queues jobs. The attachment `run.sh` is the executable run by Condor.

## Verbatim Attachment: generateArgs.sh

Source: [generateArgs.sh](mattermost_luke_may19_attachments/generateArgs.sh)

```bash
# A line will be written to args.txt with $pfn replaced by each line in pfns.txt, followed by each set of arguments
# Important: ensure there is a (blank) newline following the last entry in pfns.txt

> args.txt
while read -r pfn; do
  [[ -z "$pfn" ]] && continue
  cat >> args.txt <<EOF
$pfn boo bah blah
$pfn fah fee foe
EOF
done < pfns.txt
```

## Verbatim Attachment: pfns.txt

Source: [pfns.txt](mattermost_luke_may19_attachments/pfns.txt)

```text
root://dcgftp.usatlas.bnl.gov:1094//pnfs/blah1.root
root://dcgftp.usatlas.bnl.gov:1094//pnfs/blah2.root
root://dcgftp.usatlas.bnl.gov:1094//pnfs/blah3.root

```

## Verbatim Attachment: condor.sub

Source: [condor.sub](mattermost_luke_may19_attachments/condor.sub)

```condor
universe              = vanilla
executable            = run.sh

# Arguments for each job (`Process` is the job number, used to create unique output and error files for each job)
arguments             = $(PFNS) $(ARG1) $(ARG2) $(ARG3) $(Process)

should_transfer_files = YES
when_to_transfer_output = ON_EXIT

# Input files that your executable needs
# In this example, myMacro.C is a macro called by `run.sh`
# The file `/tmp/x509up_u102136` is where your grid proxy should be stored after you initiate it.
    # One needs to propagate the grid proxy to the worker node in order to read the file PFN
    # (When you run `voms-proxy-init -voms atlas` you should see `Created proxy in <path/to/file>`.
    # Replace `/tmp/x509up_u102136` with `<path/to/file>` if different.)
transfer_input_files  = myMacro.C, \
                        /tmp/x509up_u102136

# Assumes you have a directory `output_files` on the same level as this script
transfer_output_files = output_files

initialdir	      = .

# Assumes you have a directory `logs` on the same level as this script;
output		      = logs/runCorr_$(CORR)_$(ARG1)_$(ARG2)_$(ARG3)_$(Process).out
error		      = logs/runCorr_$(CORR)_$(ARG1)_$(ARG2)_$(ARG3)_$(Process).err
log		          = logs/condor_$(Cluster).log

# CPU, memory, disk allocation per job
request_cpus	      = 4
request_memory	      = 2GB
request_disk	      = 10MB

# `getenv=True` passes all of your shell variables to the worker node
getenv		      = True

# Grabs arguments from `args.txt` and queues the job for submission
queue PFNS,ARG1,ARG2,ARG3, from args.txt
```

## Verbatim Attachment: run.sh

Source: [run.sh](mattermost_luke_may19_attachments/run.sh)

```bash
#!/bin/bash
set -euo pipefail

WILD="$1"
ARG1="$1"
ARG2="$2"
ARG3="$3"
# `Process` (the job number) is propagated to `JOBNO``
    # (You might want `JOBNO` to propagate as an argument to your macro, or not)
JOBNO=$4

# Make the output file directory if it doesn't exist
mkdir -p "$PWD/output_files"

# For convenience define OUTDIR as the output file directory
OUTDIR="output_files"

# Export the user proxy
export X509_USER_PROXY="$PWD/x509up_u102136"

# Setup root on lxplus (change the setup script as necessary)
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_108_ATLAS_2/x86_64-el9-gcc14-opt/setup.sh
set -u

# Run your macro
root -b -l <<EOF
.L myMacro.C
myMacro("${WILD}", "${ARG1}", "${ARG2}", "${ARG3}")
.q
EOF
```

## Honest Local Checks

These checks are from comparing the downloaded attachments locally, not from Luke's chat.

- `condor.sub` passes five command-line arguments to `run.sh`: `$(PFNS) $(ARG1) $(ARG2) $(ARG3) $(Process)`.
- The downloaded `run.sh` reads `$1`, `$2`, `$3`, and `$4`; it does not read `$5`.
- The downloaded `run.sh` sets both `WILD` and `ARG1` from `$1`.
- Because this may be intentional or may be a template bug, I have not changed the attachment content. Verify the argument mapping before using these files for production jobs.

## What Was Not Verified

- I did not validate the Rucio command against a real DID.
- I did not run Condor.
- I did not verify the CERN/BNL storage instructions against current external documentation.
- I did not find any attached CC/PC stream skim list with muon trigger decisions in the May 19 Mattermost attachments.
