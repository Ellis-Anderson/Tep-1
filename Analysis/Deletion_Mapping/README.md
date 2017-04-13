###Directory for STAR deletion mapping.
####Subdirectories
#####Reads
The Reads subdirectory contains all of the raw reads to be used in STAR mapping.
Leslie's trimmed and barcode-split raw files can be found under the Day_1 Day_2 Day_3 and Day_4 directories in cyverse datashare
These directories correspond to the different days of tissue collection.
Within each Day_* directory Day*_R1.trim.paired.fq.gz contains forward split paired reads and Day*_R2.trim.paired.fq.gz contains reverse split paired reads.
Day*_all.trim.unpaired contains both R1 and R2 unpaired reads. 
#####Genome
Contains SL2.50 fasta file, seperated by chromosomes, for use in STAR mapping
#####Annotations
Contains ITAG 3.1 GFF for use in STAR mapping

