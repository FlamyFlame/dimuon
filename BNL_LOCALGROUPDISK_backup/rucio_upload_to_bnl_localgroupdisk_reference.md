# Rucio Upload To BNL LOCALGROUPDISK Reference

Source PDF: `rucio_upload_to_bnl_localgroupdisk.pdf`

Title page:

```text
Uploading datasets to the BNL LOCALGROUPDISK
Luke Mozarsky
October 21, 2025
```

This is a faithful reference-style extraction of the PDF content. It preserves the non-repeated information from the slides and notes where the original PDF repeats the same step slide while progressively adding UI actions. It does not replace `rucio_upload_to_bnl_localgroupdisk_steps.md`.

## Disk Storage For US ATLAS Users

- BNL has allocated each US ATLAS user 50 TB of permanent storage space on disk.
- The Rucio storage element (RSE) is written on the disk-storage slide as `BNL_OSG2-LOCALGROUPDISK`.
- This storage provides a way to store datasets when they are not in immediate use and download them when needed.
- Users cannot write directly to `BNL_OSG2-LOCALGROUPDISK`; instead, they must add a Rucio rule to replicate datasets from scratch disk, which is temporary disk storage.

## Recent Skims Already On Grid Scratchdisks

- Skims are stored on various scratchdisks on the grid for a specific lifetime.
- The lifetime varies from disk to disk.
- Skims can be stored in multiple datasets.
- Those datasets are grouped into a container, described in the PDF as the thing ending with `MYSTREAM`.
- Even within a single skim, datasets can be stored on different scratchdisks.

To check whether any files in the container have expired:

```bash
rucio list-files <container_name>
```

The PDF says the easiest way to add a Rucio rule is to use the Rucio UI client:

```text
https://rucio-ui.cern.ch/
```

## Rucio UI Rule Creation Steps

The PDF presents these steps cumulatively across several repeated slides titled `Uploading a recent skim to disk`.

1. Use the Rucio UI client at `https://rucio-ui.cern.ch/`.
2. Enter the scope and container name separated by a colon:

```text
user.<username>:<container_name>
```

3. Select the container and continue.

4. Choose the RSE. The upload-step slides split the name across a line break, but the RSE appears as:

```text
BNL-OSG2_LOCALGROUPDISK
```

5. Since the upload is only to one RSE, group all datasets in the container together.

6. Set a lifetime in days. Luke's slide says:

```text
I usually leave this empty
```

7. After submitting the request, navigate to your rules.

8. The rule state will say `OK` when replication to disk is complete.

## Command-Line Rule Alternative

The PDF says a Rucio rule can also be added from the command line:

```bash
rucio add-rule <scope>:<did> <number_of_copies> <RSE destination(s)>
```

The PDF defines `DID` as the name of a file, dataset, or container.

For optional arguments:

```bash
rucio add-rule -h
```

## Local Data No Longer Available On Grid

Use this section if you have data files stored locally that have expired from the grid, or files you created yourself.

First upload the files to a scratchdisk. The PDF says each user has 50 TB on:

```text
BNL_OSG2-SCRATCHDISK
```

The PDF also says a number of other scratchdisks are available. To view the list and your quota:

```bash
rucio list-account-usage <username>
```

Upload files to the scratchdisk via `rucio upload`:

```bash
rucio upload --rse <RSE> --scope <scope> --register-after-upload --recursive <local_path_to_files>
```

If you receive an error such as `Database error` when uploading a file, the PDF says Rucio likely still has some metadata or trace of the filename, and duplicate DIDs are not allowed.

Luke recommends adding the suffix `COPY` to the end of each file. Assuming all files end with `.root`, in the directory containing the files run:

```bash
for FILENAME in *; do mv $FILENAME ${FILENAME%.root}_COPY.root; done;
```

Then repeat the upload.

Make a new dataset to store these files:

```bash
rucio add-dataset <scope:new_dataset_name>
```

Attach the files to the dataset. In the directory containing the files, run:

```bash
rucio attach <scope:new_dataset_name> *.root
```

Finally, add a Rucio rule to upload this dataset to:

```text
BNL-OSG2_LOCALGROUPDISK
```

The PDF says to do this by repeating the steps from the previous several slides.

## Extraction Notes

- The PDF has 12 extracted pages, although the local file metadata reports 8 pages. The extracted pages include repeated/cumulative versions of the same upload-step slide.
- The repeated upload-step slides were consolidated above so each unique piece of information appears once.
- Plain text extraction rendered `rucio list-files` with a ligature artifact as `rucio list-ﬁles`; this Markdown uses normal ASCII `files`.
- Plain text extraction rendered `<local_path_to_files>` with a ligature artifact as `<local_path_to_ﬁles>`; this Markdown uses normal ASCII `files`.
- Plain text extraction inserted a space in the final-slide rename command as `${FILENAME%.root}_COPY .root`; I visually checked the rendered slide and used `${FILENAME%.root}_COPY.root`.
- I did not run or externally verify any command in this reference.
