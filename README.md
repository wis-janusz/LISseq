# LISseq
Small utility for analyzing LISseq data

Installation:
1. Install bowtie2 and add it to your system PATH. 
2. Clone this repository, create a Python virtual environment and install required packages using `pip install -r requirements.txt`.
3. Download a bowtie2 indexed reference genome and extract the .bt2 files to a directory.
   Default is Human / GRCh38 no-alt analysis set from https://benlangmead.github.io/aws-indexes/bowtie ,
   extracted to GRCh38_noalt_as/ in the main LISseq directory.

Usage:
1. The raw reads must be untarred to one input directory and its subdirectories but there is no need do decompress them.
2. Activate the virtual environment that you created for LISseq.
3. Execute `python LISseq.py {input_directory} {output_directory}`.
4. Identified integration sites are listed in {output_directory}/integration_sites.csv.