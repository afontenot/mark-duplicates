# mark-duplicates
Mark PCR duplicates and sequencing platform / optical duplicates in
single-end SAM files

# Introduction

[Picard](https://github.com/broadinstitute/picard) has a tool called
MarkDuplicates for marking duplicates in SAM files, but it [can only
mark optical duplicates for paired-end reads](https://gatkforums.broad
institute.org/gatk/discussion/comment/30092/#Comment_30092). This 
script is written in Python and is a very fast program for doing the 
same for single-end reads.

# Usage

[PyPy](https://pypy.org/) is strongly recommended for its performance
improvement over standard Python. It's several orders of magnitude 
faster for this program.

The input and output of the program are uncompressed SAM files in 
fastq coordinates. The program could easily be extended to add support
for compression and decompression (though full BAM support would be 
more complicated), but if you have space concerns I suggest 
compression on a file-system basis.

To run the program:

    pypy mark_duplicates.py input.sam output.sam pixels

Note that `pixels` is a required argument. It is a radius (in **qseq**
pixel coordinates) within which reads will be marked as optical 
duplicates.

Similar to the Picard tool, this program marks duplicates with the
following notation:

 * `DT:Z:SQ` for optical duplicates
 * `DT:Z:LB` for PCR duplicates

The duplicate flag is also set. Note that any candidate optical
duplicate is also detected as a PCR duplicate but only marked as an
optical duplicate. If your analysis requires removing one but not the 
other you should take this into account.

If you want either of these duplicates removed instead, this can
be accomplished by changing only one line of the program. (Readability
and making the program easy to change were major goals.)

# Assumptions

This program was written to be used for a specific biology paper, and
it makes the following assumptions, which may or may not be safe for
your project.

 * Input is a single-end, uncompressed SAM file in fastq coordinates.
 * Input is sorted by mapping position.
 * Each group of reads (reads with the same mapping position) is small
   enough to fit into memory. (This vastly speeds up the program.)
 * In the output, if `n` reads are detected to be duplicates, only the
   `n-1` lower [sequence quality](http://www.drive5.com/usearch/manual/
   quality_score.html) reads will be marked as duplicates. 
   The best read will remain unmarked.
 * In the output, if a gene is detected to be a duplicate, the 
   corresponding duplicate flag should be set.
 * The reads in the output file may be resorted.

# License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
