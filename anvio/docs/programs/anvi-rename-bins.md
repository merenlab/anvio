This program **creates a new %(collection)s from the %(bin)ss in another collection with specific guidelines.** This is especially helpful when you want to merge multiple collections later or share your project with someone, and you want all of your bins to have nicer names than the default `bin_01`, `bin_02`, etc. based on the order you binned them in.

So let's take a look at what this program can do with a simple example.

### Example 1: Renaming all bins in a collection

Let's say you have a collection called `MY_COLLECTION`, which has four bins: `really`, `bad`, `bin`, and `names`. These names just won't do, so let's get to renaming. To rename all of my bins and put them into a collection called `SURFACE_OCEAN_SAMPLES`, you could run

{{ codestart }}
anvi-rename-bins -c %(contigs-db)s \
                 -p %(profile-db)s \
                 --prefix SURFACE_OCEAN \
                 --collection-to-read MY_COLLECTION \
                 --collection-to-write SURFACE_OCEAN_SAMPLES \
                 --report-file rename.txt
{{ codestop }}

And voila! Now you have a second collection named `SURFACE_OCEAN_SAMPLES` that contains your four bins, now named  `SURFACE_OCEAN_Bin_00001`, `SURFACE_OCEAN_Bin_00002`, `SURFACE_OCEAN_Bin_00003`, and `SURFACE_OCEAN_Bin_00004`. The order that the numbers are in represents the quality of the bin as a MAG, given by the completion minus redunancy.

The file `rename.txt` is just a tab-delimited file that contains a summary of your renaming process. The first column has the original name of the bins that you renamed, the second has their new names, and the remaining columns contain information about those bins (like their completion, redundency, and size).

### Example 2: Separating out the MAGs

Okay, but what if you want to label your MAGs separately from your bins? You don't like `SURFACE_OCEAN_bin_00004` since it only has a completition stat of 50 percent, and you're not sure if you want to include `SURFACE_OCEAN_bin_00003`  since it has 50 percent redundency. How can you differenciate these iffy bins in your collection?

Here is the solution:

{{ codestart }}
anvi-rename-bins -c %(contigs-db)s \
                 -p %(profile-db)s \
                 --prefix SURFACE_OCEAN \
                 --collection-to-read MY_COLLECTION \
                 --collection-to-write SURFACE_OCEAN_MAGS \
                 --report-file rename.txt \
                 --call-MAGs \
                 --min-completion-for-MAG 70
{{ codestop }}

Now, the collection `SURFACE_OCEAN_MAGS` will include  `SURFACE_OCEAN_MAG_00001`, `SURFACE_OCEAN_MAG_00002`, `SURFACE_OCEAN_MAG_00003`, and `SURFACE_OCEAN_Bin_00004`. These are exactly the same bins that the collection contained before, but now the names differenciate the wheat from the chaff.

Now, let's make that same collection (still called `SURFACE_OCEAN_MAGS`) that doesn't include `SURFACE_OCEAN_Bin_00003` as a MAG, since the redundency is too high for what we want to look at right now.

{{ codestart }}
anvi-rename-bins -c %(contigs-db)s \
                 -p %(profile-db)s \
                 --prefix SURFACE_OCEAN \
                 --collection-to-read MY_COLLECTION \
                 --collection-to-write SURFACE_OCEAN_MAGS \
                 --report-file rename.txt \
                 --min-completion-for-MAG 70 \
                 --max-redundancy-for-MAG 30 \
                 --call-MAGs
{{ codestop }}

Now `SURFACE_OCEAN_MAGS`   will include  `SURFACE_OCEAN_MAG_00001`  `SURFACE_OCEAN_MAG_00002`,  `SURFACE_OCEAN_Bin_00003`, and `SURFACE_OCEAN_Bin_00004`.

You also have the option to only classify bins above a certain minimum size as MAGs.

### Example 3: An example use case in a workflow

For an example use case, on [this page](http://merenlab.org/tutorials/infant-gut/#renaming-bins-in-your-collection-from-chaos-to-order), anvi-rename-bins is used to create a new collection called `MAGs` that contains differenciates bins that have a completion stat of more than 70 percent, and renames all of those bins with the prefix `IGD` (which stands for infant gut dataset).
