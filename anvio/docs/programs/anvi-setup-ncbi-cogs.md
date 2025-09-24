This program **downloads and organizes a local copy of the data from NCBI's [COGs database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/) for use in function annotation.** This program generates a %(cogs-data)s artifact, which is required to run the program %(anvi-run-ncbi-cogs)s.

### Set up COGs data

{{ codestart }}
anvi-setup-ncbi-cogs
{{ codestop }}

{:.warning}
If possible, we recommend you multithread this program with the `--num-threads` parameter to make it faster.

If you already have a %(cogs-data)s artifact and are trying to redownload this data, run

{{ codestart }}
anvi-setup-ncbi-cogs --reset
{{ codestop }}

### Choosing a different database version

{{ codestart }}
anvi-setup-ncbi-cogs --cog-version COG14
{{ codestop }}

Not sure which versions of %(cogs-data)s are available? You can type something random after the `--cog-version` parameter to see the options.

### Manually downloading the raw data from the NCBI

Some users may need to download the raw data from the NCBI manually, instead of letting `anvi-setup-ncbi-cogs` to do it for them. There may be multiple reasons for this. One reason is the NCBI servesr. Unfortunately we learned the hard way that the NCBI servers are not always reliable, and sometimes the download process gets interrupted, leading to checksum errors, especially for users who have exceptionally slow internet connections. Another reason is that some users may be working in a computing environment that does not have direct access to the internet, and therefore cannot download the data directly from the NCBI. In such cases, users can download the data on a different computer that has (fast) internet access, and then transfer the raw data files to the computer where they want to set up the COGs database.

One issue here is to know exactly what to download, and where to put those files depending on the `--cog-data-dir` or `--cog-version` parameters. Luckily, the `anvi-setup-ncbi-cogs` program can help with that if you were to run it with the `--dry-run` flag on the computer where you wish to setup the COGs database. Here are some examples and the output messages the program will generate for you (I am running these commands on my own comptuer that uses `anvio-dev`, but it will generate specific instructions for your own comptuer):

```
anvi-setup-ncbi-cogs --dry-run

COG version ..................................: COG24
COG data source ..............................: The anvi'o default.
COG base directory ...........................: /Users/meren/github/anvio/anvio/data/misc/COG

DRY RUN MODE WILL NOW TELL YOU THINGS AND QUIT
===============================================
The following information will tell you what anvi'o would have done to acquire
the raw files from the NCBI to set them up on your computer. Using this
information you should be able to download these files, and move them manually
to the locations shown below. Please remember that the data shown below will
depend on teh COG version of your choosing (or, if you haven't specified any,
the default COG version). Before you start following the instructions below,
make sure the following directory exists:

    - /Users/meren/github/anvio/anvio/data/misc/COG/COG24/RAW_DATA_FROM_NCBI

Then execute the following steps like a robot 🤖

Download this file ...........................: https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.cog.csv
 -> and move it here .........................: /Users/meren/github/anvio/anvio/data/misc/COG/COG24/RAW_DATA_FROM_NCBI/cog-24.cog.csv

Download this file ...........................: https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab
 -> and move it here .........................: /Users/meren/github/anvio/anvio/data/misc/COG/COG24/RAW_DATA_FROM_NCBI/cog-24.def.tab

Download this file ...........................: https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab
 -> and move it here .........................: /Users/meren/github/anvio/anvio/data/misc/COG/COG24/RAW_DATA_FROM_NCBI/cog-24.fun.tab

Download this file ...........................: https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/checksums.md5
 -> and move it here .........................: /Users/meren/github/anvio/anvio/data/misc/COG/COG24/RAW_DATA_FROM_NCBI/checksum.md5.txt

Download this file ...........................: https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/COGorg24.faa.gz
 -> and move it here .........................: /Users/meren/github/anvio/anvio/data/misc/COG/COG24/RAW_DATA_FROM_NCBI/COGorg24.faa.gz

Once you have all these files in those locations, `anvi-setup-ncbi-cogs` should
not attempt to download anything, but process the existing files at these
specific locations to setup NCBI COGs on your system in peace.
```

Or if I were to use a specific version and a different COG data directory, I would have gotten different instructions:

```
anvi-setup-ncbi-cogs --dry-run --cog-data-dir /tmp/test --cog-version COG20

COG version ..................................: COG20
COG data source ..............................: The command line parameter.
COG base directory ...........................: /tmp/test

DRY RUN MODE WILL NOW TELL YOU THINGS AND QUIT
===============================================
The following information will tell you what anvi'o would have done to acquire
the raw files from the NCBI to set them up on your computer. Using this
information you should be able to download these files, and move them manually
to the locations shown below. Please remember that the data shown below will
depend on teh COG version of your choosing (or, if you haven't specified any,
the default COG version). Before you start following the instructions below,
make sure the following directory exists:

    - /tmp/test/COG20/RAW_DATA_FROM_NCBI

Then execute the following steps like a robot 🤖

Download this file ...........................: ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv
 -> and move it here .........................: /tmp/test/COG20/RAW_DATA_FROM_NCBI/cog-20.cog.csv

Download this file ...........................: ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
 -> and move it here .........................: /tmp/test/COG20/RAW_DATA_FROM_NCBI/cog-20.def.tab

Download this file ...........................: ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/fun-20.tab
 -> and move it here .........................: /tmp/test/COG20/RAW_DATA_FROM_NCBI/fun-20.tab

Download this file ...........................: ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/checksums.md5.txt
 -> and move it here .........................: /tmp/test/COG20/RAW_DATA_FROM_NCBI/checksum.md5.txt

Download this file ...........................: ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.fa.gz
 -> and move it here .........................: /tmp/test/COG20/RAW_DATA_FROM_NCBI/cog-20.fa.gz

Once you have all these files in those locations, `anvi-setup-ncbi-cogs` should
not attempt to download anything, but process the existing files at these
specific locations to setup NCBI COGs on your system in peace.
```

### Always getting checksum errors?

Sometimes the NCBI servers drop incoming connections, leading to incomplete file downloads and errors like the following:

```
Config Error: Something went wrong with your download :/ The checksum we calculated for
              `cog-20.cog.csv` anvi'o just finished downloading does not match to the checksum
              provided by the NCBI. This is most likely due to an interrupted download, as the
              NCBI servers often prematurely end data transfers. Please try running the same
              command again with the `--reset` flag.
```

If you have tried re-running `anvi-setup-ncbi-cogs` multiple times but are always getting checksum errors and are about to lose your mind, here is a set of commands that you can follow to manually download the data for the **2020 release of COGs** without having to go through the setup program every time.

First, you will need to move to the directory where anvi'o expects to find the COG files. This location will depend on where conda and anvi'o are installed on your computer, but if you have the anvi'o environment loaded in your terminal, you can easily get there by running the following:

```
cd $CONDA_PREFIX/lib/python3.10/site-packages/anvio/data/misc/COG/COG20/RAW_DATA_FROM_NCBI/
```

{:.warning}
This example is for the 2020 release of COGs, but you can figure out which directory to go based on the COG version on your computer by running `anvi-setup-ncbi-cogs --dry-run` as described above.

The files that anvi'o needs to see in that folder are the following:
```
checksum.md5.txt cog-20.def.tab   fun-20.tab
cog-20.cog.csv   cog-20.fa.gz
```

Since you have already tried running `anvi-setup-ncbi-cogs` so many times, probably there are some of those files already in there. But the checksums of those files need to match those that are listed in the `checksum.md5.txt` file. For instance, if you look for `cog-20.cog.csv` inside the checksum file:

```
grep cog-20.cog.csv checksum.md5.txt
```

You will see the following line:

```
1bed944a61e0ec404669361fb69ae52d  cog-20.cog.csv
```

which indicates that the file's checksum should match exactly to `1bed944a61e0ec404669361fb69ae52d`. If you run `md5sum cog-20.cog.csv`, you should see that exact string. If you don't see the same thing, it means the file has been incompletely downloaded, so it needs to be downloaded again. You can do it like this:

```
rm -rf cog-20.cog.csv

curl -O https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.cog.csv

md5sum cog-20.cog.csv
```

If you don't want to compare the strings manually, you can use the `diff` program to check that they are the same, like this:

```
grep cog-20.cog.csv checksum.md5.txt > expected_checksum
md5sum cog-20.cog.csv > observed_checksum
diff expected_checksum observed_checksum
```

If the two checksums match, you won't see any output from the `diff` command.

Once you get a copy of the file with an exactly matching MD5 checksum, you can move on.

You should run `md5sum` on every file listed above (except for `checksum.md5.txt`), and check if it matches the corresponding string inside `checksum.md5.txt`. For any file with a non-matching MD5 checksum, you should download it using `curl` as we did above:

```
rm -rf [FILENAME THAT DOES NOT MATCH]
curl -O https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/[FILENAME THAT DOES NOT MATCH]
```

(make sure you change the file name at the end of the path to match the file that you need)

After you have all the files with matching checksums, you can leave the data folder, and then re-run `anvi-setup-ncbi-cogs`, which should now work perfectly using the manually downloaded files:

```
cd
anvi-setup-ncbi-cogs
```
