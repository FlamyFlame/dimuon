# Rucio Upload Steps To BNL LOCALGROUPDISK

Source PDF: `rucio_upload_to_bnl_localgroupdisk.pdf`, titled "Uploading datasets to the BNL LOCALGROUPDISK" by Luke Mozarsky, October 21, 2025.

This file extracts only the upload/replication steps from the PDF. I did not verify these commands against current Rucio documentation or run them.

## Context From The PDF

- BNL provides US ATLAS users permanent disk storage, described in the PDF as 50 TB per user.
- The PDF says users cannot write directly to LOCALGROUPDISK.
- The workflow is to replicate datasets to LOCALGROUPDISK by adding a Rucio rule.
- The PDF refers to BNL LOCALGROUPDISK with slightly different spellings in different places. The step slides and final local-upload slide use `BNL-OSG2_LOCALGROUPDISK`; one context slide text extraction shows `BNL_OSG2-LOCALGROUPDISK`.

## Workflow A: Recent Skim Already On The Grid

Use this when the skim still exists on grid scratchdisk storage.

1. Check whether any files in the container have expired:

```bash
rucio list-files <container_name>
```

2. Open the Rucio UI client:

```text
https://rucio-ui.cern.ch/
```

3. Start adding a Rucio rule in the UI.

4. Enter the scope and container name separated by a colon:

```text
user.<username>:<container_name>
```

5. Select the container and continue.

6. Choose the RSE:

```text
BNL-OSG2_LOCALGROUPDISK
```

7. Since the upload is to one RSE, group all datasets in the container together.

8. Set a lifetime in days. Luke's slide says he usually leaves this empty.

9. Submit the request.

10. Navigate to your Rucio rules. The state will say `OK` when replication to disk is complete.

## Command-Line Alternative For Adding A Rule

The PDF says a rule can also be added from the command line:

```bash
rucio add-rule <scope>:<did> <number_of_copies> <RSE destination(s)>
```

Notes from the PDF:

- `DID` is the name of a file, dataset, or container.
- Run this for optional arguments:

```bash
rucio add-rule -h
```

## Workflow B: Local Data No Longer Available On The Grid

Use this if files are local because they expired from the grid, or because you created the files yourself.

1. Choose a scratchdisk. The PDF says each user has 50 TB on `BNL_OSG2-SCRATCHDISK`, with other scratchdisks also available.

2. View available scratchdisks and quota:

```bash
rucio list-account-usage <username>
```

3. Upload local files to scratchdisk:

```bash
rucio upload --rse <RSE> --scope <scope> --register-after-upload --recursive <local_path_to_files>
```

4. If upload fails with an error such as `Database error`, the PDF says Rucio may still have metadata or a trace of the filename, and duplicate DIDs are not allowed.

5. Luke recommends adding the suffix `COPY` to each file. For files ending in `.root`, run this in the directory containing the files:

```bash
for FILENAME in *; do mv $FILENAME ${FILENAME%.root}_COPY.root; done;
```

6. Repeat the upload command after renaming.

7. Make a new dataset to store these files:

```bash
rucio add-dataset <scope:new_dataset_name>
```

8. Attach the uploaded files to the dataset. In the directory containing the files, run:

```bash
rucio attach <scope:new_dataset_name> *.root
```

9. Add a Rucio rule to replicate this dataset to `BNL-OSG2_LOCALGROUPDISK` by following Workflow A.

## Visual Extraction Note

The final-slide rename command was checked visually because plain text extraction inserted a space before `.root`. The visually verified command is:

```bash
for FILENAME in *; do mv $FILENAME ${FILENAME%.root}_COPY.root; done;
```
