Reports sequences stored in a %(bam-file)s file that cover one of more specific nucleotide positions in a reference.

### Basic mode of operation

Assume you wish to recover reads stored in one or more BAM files, where the matching reads contain at least one nucleotide position that align to a nucleotide position `nucleotide_position_N` in a contig `contig_name_X`. In that case, the user would first generate a two column TAB-delmited file, for example `positions_for_linkmers.txt` with no header line,


<table>
  <tbody>
    <tr>
      <td> contig_name_X </td>
      <td> nucleotide_position_N </td>
    </tr>
  </tbody>
</table>

And run the program this way to recover the short reads from this way:

```
anvi-report-linkmers --contigs-and-positions positions_for_linkmers.txt \
                     -i SAMPLE_01.bam SAMPLE_02.bam SAMPLE_03.bam (...) \
                     -o linkmers.txt
```

The user can define multiple contigs in the input file, and one or more nucleotide positions for each one of them:

<table>
<tbody>
<tr>
      <td> contig_name_X </td>
      <td> nucleotide_position_01,nucleotide_position_02,nucleotide_position_03</td>
</tr>
<tr>
      <td> contig_name_Y </td>
      <td> nucleotide_position_04 </td>
</tr>
<tr>
      <td> contig_name_Z </td>
      <td> nucleotide_position_05,nucleotide_position_06 </td>
</tr>
</tbody>
</table>

The resulting %(linkmers-txt)s would include all short reads that match any of these critera

### Complete or incomplete links?

Using the `--only-complete-links` flag, the user can enforce whether only complete links should be reported where each reported short read must cover each nucleotide position for a given contig.

Please note that if the nucleotide positions chosen for a given contig are too distant from each other given the short read length, zero reads may satisfy the complete links criterion.

Having complete links, however, will enable [oligotyping](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12114) analyses on **metagenomic reads** through the anvi'o program %(anvi-oligotype-linkmers)s.

### See this program in action

[http://merenlab.org/2015/12/09/musings-over-commamox/](http://merenlab.org/2015/12/09/musings-over-commamox/)