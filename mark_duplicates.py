# mark_duplicates.py: mark duplicates in single-end SAM files
# Copyright (C) 2019  Adam Fontenot
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

from __future__ import print_function
import os
import sys


# quality for a base = the ascii ordinal - 64
# we define quality for a read to be the sum of base qualities
def quality(qualitystring):
    totalquality = 0
    for i in qualitystring:
        totalquality += ord(i) - 64
    return totalquality


# get square of distance between points for fast compare
def euclid_sq(x, y):
    return (int(x[4]) - int(y[4])) ** 2 + (int(x[5]) - int(y[5])) ** 2


# convert qseq pixel coordinates to fastq format (squared)
def qseq_to_fastq(qseq_pixels):
    return (10 * qseq_pixels + 1000) ** 2


# get sequence for key based sorting of gene groups
def sequence(read):
    return read[9]


# get number of tile for key based sorting of gene groups
def tile_num(read):
    return int(read[0].split(':')[3])


def main():
    if len(sys.argv) != 4:
        print("syntax: {} input.sam output.sam pixels".format(sys.argv[0]))
        return
    f = open(sys.argv[1], "r")
    fsize = os.path.getsize(sys.argv[1])
    o = open(sys.argv[2], "w")
    maxdistance = qseq_to_fastq(int(sys.argv[3]))
    
    nextprint, linenumber = 0, 0
    groupline = ""
    while True:
        newline = ""
        if groupline == "":
            newline = f.readline()
            if newline == "":
                break
        else:
            newline = groupline
        linenumber += 1
        if linenumber > nextprint:
            nextprint = linenumber + 10000
            print('\r' + str(100 * f.tell() // fsize) + '% complete', end='')
            sys.stdout.flush()
        read = newline.split('\t')
        # skip bad reads
        if read[2] == '*':
            continue
        readgroup = [read]
        
        ## loop 1: collect a group of reads
        while True:
            pos = f.tell()
            groupline = f.readline()
            if groupline == "":
                break
            groupread = groupline.split('\t')
            # group ends when we leave gene or location
            # groupline is stored to use next time through the outer loop
            if (read[2] != groupread[2] or read[3] != groupread[3]):
                break
            linenumber += 1
            readgroup += [groupread]
            
        # no need to check for duplicates if there's only one read in group
        if len(readgroup) == 1:
            o.write('\t'.join(readgroup[0]))
            continue
        if len(readgroup) > 1000:
            print("\nbeginning large group of size", len(readgroup))
        
        isdupe = [False] * len(readgroup)
        isopticaldupe = [False] * len(readgroup)
        
        # presort by sequence and tile so we can break out quickly later
        readgroup.sort(key=lambda x: (sequence(x), tile_num(x)))
        # cache an array for read quality
        groupquality = [quality(i[10]) for i in readgroup]
        
        ## loop 2: detect ordinary duplicates
        for i in range(len(readgroup)):
            # if i is a dupe, there's a better quality non-dupe, so skip
            if isdupe[i]:
                continue
            # only mark the lower quality read as a duplicate
            for j in range(i+1,len(readgroup)):
                # don't mark reads with different sequences as duplicates
                if readgroup[i][9] != readgroup[j][9]:
                    break
                if groupquality[i] >= groupquality[j]:
                    isdupe[j] = True
                    readgroup[j][1] = str(int(readgroup[j][1]) | 0x400)
                else:
                    isdupe[i] = True
                    readgroup[i][1] = str(int(readgroup[i][1]) | 0x400)
        
        # cache readnames for finding opticals
        groupreadnames = [i[0].split(':') for i in readgroup]
        
        ## loop 3: find optical duplicates
        for i in range(len(readgroup)):
            for j in range(i+1,len(readgroup)):
                # don't mark reads with different sequences as duplicates
                if readgroup[i][9] != readgroup[j][9]:
                    break
                # if tile number doesn't match, skip optical checking
                if groupreadnames[i][3] != groupreadnames[j][3]:
                    break
                # if lane number doesn't match, not an optical dupe
                if groupreadnames[i][2] != groupreadnames[j][2]:
                    continue
                distance = euclid_sq(groupreadnames[i], groupreadnames[j])
                if distance < maxdistance:
                    # only mark the lower quality read as a duplicate
                    if groupquality[i] >= groupquality[j]:
                        isopticaldupe[j] = True
                    else:
                        isopticaldupe[i] = True
                        
        ## annotate output
        for i in range(len(readgroup)):
            if isopticaldupe[i]:
                o.write('\t'.join(readgroup[i])[:-1] + " DT:Z:SQ\n")
            elif isdupe[i]:
                o.write('\t'.join(readgroup[i])[:-1] + " DT:Z:LB\n")
            else:
                o.write('\t'.join(readgroup[i]))

    print("\r100% complete")

if __name__ == "__main__":
    main()
