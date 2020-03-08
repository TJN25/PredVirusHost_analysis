# PredVirusHost
Predicts the host domain of contigs from viromes by comparing proteins from a number of 
contigs to hmm models and then scores each contig based on similarity to those models.



Dependancies:
- HMMER (http://hmmer.org/download.html)
- R
- dplyr (R package)



The files that will be found are:
-predvirushost.sh
	This is the wrapper used to carry out the prediction of the viral hosts. 
-host_scoring.R
	Takes the output of the hmmsearches and scores each contig. It is called by the 
	predvirushost.sh script.	
-*.hmm
	There are 3 profile-HMM files (arVOG, baPOG and euVOG) and these are used in the
	hmmsearch step
_*model_scores.txt
	For each of the sets of profile-HMMs (arVOG, baPOG and euVOG) there are scores that are used to weight each 
	model depending on how well the model can discriminate between the true hosts and 
	the other hosts. These scores range from 0-1 and are used to adjust the bit scores 
	in the host_scoring.R script.

Options:
The input requires that a fasta file (protein) is specified and that the format of the 
protein headers is also specified. This is required for files with more than one contig as
it is necessary to be able to idenitfy the protein name and contig name so that the results
of the hmmsearch can be matched to the correct contig for scoring. There are default options
for this (MGRAST, PROKKA, Refseq, Other). If 'Other' is selected then the format of the 
header will need to be specified by giving the character that separates the contig name from
the rest of the protein name along with the position of this character in the header.

Example 1.

>80818167_1_425_+
XQXFETAARFGSIDDIRYKRPHHTRIVKDEILGTIRQVLTITPQEITLSRGDIDIIVRQAKVCQVTKALKGIDTIIAYDRGFEGEGRSHWRGVDIVPLHSDHVDKLLRKIVDFDTLANVLECLRNTCREISTHNVVGAT
>80818167_525_896_-
VAKTETKTTKIPDDVKQYVDIIIDALSKALETGNFTIELTYDNIKIKTGDINAYIDYKKLYIMYKDIDIIFSDYVAMVKVYNDYYSDPQIYYAHEYWAILEELHATAIDKVKERIKKALTKI
>80818167_942_1050_-
ALCSCVNMALQNANTLADMKYVPHDVLDTDKGHG

The contig name is '80818167'
This is separated by a '_' (-d _)
This is from the start (-s 0)
This ends 3 '_' from the end of the header (-e 3) 

Example 2.

>ABP73391.1_hypothetical_protein_[Acidianus_bottle-shaped_virus_complete_genome.]
MTTTHDTNTKKLKYQFHTIHSQRIMTTVTQKPFTASPYIFSTTLRTTQTDGNNAINSHSHTQAGYNNSSERFLYLICTYIT

>ABP73392.1_hypothetical_protein_[Acidianus_bottle-shaped_virus_complete_genome.]
MNTQNIVEKTHEKFSSTWSYLEQTSTDLKQNLRSRLLGLDSLLRELYKQETDEKKKMYIFYALMYVDSALDNLDYMNPEHDPFAVFQAKWAVQKAFEIFSDLFKDYFKSD

>ABP73393.1_hypothetical_protein_[Acidianus_bottle-shaped_virus_complete_genome.]
MESLIKAIREEFISIFSLLKKPHKFTLEELRIRLRKLHNLLLQLFELEYDINRKVEVKYALEFVDSVYNQIDYVIPETLIPFAKEQLKYAVYRFNLYYT

The contig name is 'Acidianus_bottle-shaped_virus_complete_genome'
This is separated by a '[' (-d [)
This is at the first occurrence of the '[' (-s 1)
This continues to the end of the header (-e 0)

Formatting the input file:
If the contig labels in the output table are not correct then use the following for 
the input file.

##FILE-NAME.faa
>protein_id[contig_name]
protein_sequence

>protein_id[contig_name]
protein_sequence

>protein_id[contig_name]
protein_sequence


##
Make sure that any spaces have been removed or replaced. 
Use:
 -f Other -d [ -s 1 -e 0 
along with any other options.




The help page shows:


PredVirusHost.sh: compares proteins from a number of contigs to hmm models and then scores
each contig based on similarity to those models.
Version 3.0 2018
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
usage: [options] <filename> <file type>
 
Basic options:
  -i : <filename> the protein fasta file with protein identifiers not containing spaces
  -o : <output file> the output file name (default is the input file name)
  -f : <file type> the format of the fasta file (MGRAST, PROKKA, Refseq, Other)
  -c : <Number of cores> The number of cores to use in the hmmsearch step
  -p : <Number of proteins> contigs with this number of proteins or more will be analysed with the pipeline 
Options if "Other" file type is selected:
  -d : <delim> is the character that separates the contig name from the protein name/detail
  -s : <start_pos> is the number of delimiters that you must pass to be left with just the contig name.
  -e : <end_pos> is the number of delimiters from the end needed to be left with just the contig
  For MG-RAST files the format would be -d _ -s 0 -e 3
 
Trouble shooting:
  -k : include in order to keep the temp files to see where there is a problem


The output file will contain a tab-delimited table.

##Example table
genome	archaeal_score	phage_score	eukaryotic_score	call	difference.percentage	protein.counts
contig100010	38.493	0	0	archaeal	100	1
contig309	34.692	30.295	0	archaeal	12.67	5
contig3338	39.298	95.024	0	phage	58.64	3
contig273	182.7	44.8	584.002	eukaryotic	68.72	6

Columns are:
-Genome
	The contig name
-*_score
	archaeal, phage, eukaryotic
	The sum of the adjusted bit scores for the contig for each set of models.
-call
	The host that the contig has been assigned to.
-difference.percentage
	The percentage difference between the score of the assigned host and the next
	highest scoring host.
-protein.counts
	The number of proteins in the contig.
	




