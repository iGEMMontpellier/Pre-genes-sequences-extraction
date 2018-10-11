# Pre-genes-sequences-extraction
The goal of this code is to help  extract pre-gene sequences form a full genome in the goal to do bioinformatic treatment on those.


Killy take as input file a genbank file (.gb or gbk) standar for full annoted genome, and extract  all pre-gene sequence with a size from 25 to 100 bp (declared variables at the begining, can be modified to fit your need) and check for overlaps. If the pre-gene sequence is not in the direct sense (strand=-1) it take the reverse complement of the sequence. finally it creat a .fasta file to use for MEME suite (http://meme-suite.org/tools/meme). with each sequence being name "prom" + "N = order on the genome" + "strand" (nothing if direct, -1 if reverse complement) + "locus tag". 
As the sequence are all in the good order you don't need to do an analisis on both strands with MEME, saving you some times.

To use MEME upload the fatsta file, and enter you parameters (number of motif you want to find...), if you have a lot of sequences search for motif with a minimun of 30 sites per 1000 sequence, to avoid aberations.

Then you got the result, with a meme and a mast, geting you the localisation of each motif on your sequences. If like me you want an histogramme you can use Cibo, for this open the .xml mast file and download it (keep the name mast.xml), and run cibo in the same folder. and tadam a corect histogramme of the localisation of the sequence. 
