# LISseq
Small utility for analyzing LISseq data

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC_BY--NC--ND_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

Installation:
1. Download the latest release and unzip to a directory.
2. Execute `./setup.sh`
3. Install bowtie2 and add it to your system PATH. 
4. Download or build a bowtie2 reference genome index and extract the .bt2 files to a directory.
   Default is Human / GRCh38 no-alt analysis set from https://benlangmead.github.io/aws-indexes/bowtie ,
   extracted to subdirectory GRCh38_noalt_as/ in the directory where you placed LISseq.

Usage:
1. The raw reads must be untarred to one input directory and its subdirectories, but there is no need do decompress them.
2. Execute `./run.sh {options} {input_directory} {output_directory}`.
4. Identified integration sites are listed in {output_directory}/integration_sites.csv.
