# Project crRNA tracker
My pipeline that locates potential Crispr primer regions for eDNA analysis. 

This walkthrough shall explain how the crRNA selection pipeline works.

First we select two species: our species of interest and a close relative. In this case we have the complete mitocho\
ndrion of Atlantic Salmon and Atlantic trout.

For the record, doing a quick global alignment between the two mitochondria using Waterman-Smith gives us 83.9% simil\
arity between the two sequences. this is because both of them have extended tailing ends so doing a local alignment \
gives us 94.7% similarity.

The template for searching primers is as follows; the site that we are looking should be structured as follows:


PAM_SITE ----> CRISPR_RNA ----> CAS12

In base pair format: TTTV-NNNNNNNNNNNNNNNNNNNN

Description: The **PAM site** should be begin with three Ts and end in a non-T base pair.
             The **CRISPR RNA** (crRNA) site  (which will be the primer template) should be 20 bp long without containin\
g a PAM pattern (a TTT pattern)

Now what we would like to find out is a PAM_guide site that would be distinctive to salmon and not align with trouts\
. In the paper we used as a guide, the PAMsite they selected was:

    Salmon Site
    TTTC TACCCTCCAAAACCCCTATC

    ####  #              ##
    PAM site       Polymorphic bases

    Trout Site:
    TTCC TCCCCTCCAAAACCCCCGTC

Therefore when comparting the two sequences,we want to find a site that has:

-The correct PAM site structure.
-A guide RNA that is different enough so as, if the PAM site is misidentified, the crRNA wont bind to the other site\
.

So the short version of the steps we take to do this goes as follows.

1. Download the two files containing the complete mitochondrion sequence.
2. With the species of interest, identify all sequences that are potential crRNA sites, isolate and index them.
3. Do a local alignment between the crRNA sequences and the mitochondrial genome of the compared species.
4. Find which is the lowest level of similarity between the aligned forms (out of 20) and select those sequences.
   *might also be useful to make sure that the two sequences belong to the same gene region, just to make sure we ar\
e looking at the same site.
5. For one extra safety step, if there are more than one sites that have the same level of similarity, look at the P\
AM sites of these sequences and select for the ones that are not initiated by a PAM site.

=======
**The pipeline (at some point I will compile everything together)**

**Starting material:**
Subj_FASTA
Ref_FASTA

Extra material:
Subj_cds_FASTA ---> FASTA file that contains complete information about the coding regions

----

**clip_guide_sites.pl**
Input file: Subj_FASTA
Function: Parses the FASTA file, identifies PAM sites and creates an enumerated fasta file of potential crRNA sequen\
ces with unique identifiers. An extra file is created which gives the location of each sequence on the mitochondrion\
.
Output file: subject_keys.fa
             subject_map.txt

**watering_sequences.sh**
Input file: subject_keys.fa
            reference_FASTA
Function: Performs a local alignment between the potential crRNA sites and the reference site
Output file: reference_vs_subject.water

**parse_water_check_pams.pl**
Input file: reference_vs_subject.water OR a grepped file that only contains the lowest level of similarity (for exam\
ple in the salmon vs trout) the sequences had 17/20 bases matching so I grepped for that by doing grep -A7 -B7 17/20\
 to capture the correct values. And on that note there was one more sequence that had 16/20 level of similarity so t\
hat could also be a potential site.

**parse_water_check_pams_min.pl**
An updated version of parse_water_check_pams.pl that includes the extra functions of selecting and printing out on a separate file only the sequences that have the minimum value out of 20 (so very simple for now though I have some plans in mind) on a separate file along with a file that includes all the key values with the associate sequences for reference, a file ending in "..water_map.txt".

Function: - reads in the water file and parses the reference sequences that matched with the subject sequence, picki\
ng the same unique numerical identifies for their FASTA id.
          - aligns those sequences with the file

Output file: output_file.txt


