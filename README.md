# PredVirusHost
Predicts the host domain of contigs from viromes
The files that will be found are:
-Virus-Host_Matcher.sh
	This calls other scripts, makes some simple adjustments to the files and sorts the final results into a folder.	
-remove_small_contigs.py
	This looks at how many proteins are in each contig and removes any contigs from the fastafile if they
	contain less than the number indicated by the user (default is 5 proteins or more).
-match_proteins_to_contigs.py 
	This requires that the user input be correct. This takes the whole protein name, and separates out the part that will 
	help define which contig the protein belongs to. For MG-RAST data, the charater where the split (delim) will happen is _ and the 
	contig name is from the first part (start pos) (0) through until (end pos) the 3rd _ is found (3) (-g _ -s 0 -e 3). For Refseq files (once all 
	spaces have been removed) it would be (-g [ -s 1 -e 0).
-find_and_replace_space_with_tab.py
	The output of hmmsearch contains a changing number of spaces. This is needed so the other scripts can easily work with the files.
-generate_output_file.py
	Combines the three hmm outputs into one fiile, where the best match for each protein to each host domain is kept.
-significant_matches_output.py
	Provides a file where only proteins that have matched something are kept.
-output_with_contig_names.py
	Appends the contig/species name onto the end of the output file.
-format_for_scoring.py
	Generates a file that allows R to quickly score the contigs.
-score_matches_for_contigs.R
	Scores the contigs, and includes an overall p-value and the mean and median e-values as well as the number of 
	proteins on each contig.
-format_output_scores.py
	Formats so the output contains headers and is tab-delimited to allow further analysis to be easily done.
test.faa is a test file for checking everything is working and the Sample_Results folder contains 
the expected output from running Virus-Host_Matcher.sh.


The files containing the significant models are:
arVOG_sig.hmm
euVOG_sig.hmm
baPOG_sig.hmm

The files containing all of the models are:
arVOG-all.hmm
euVOG_all.hmm
baPOG_all.hmm


There are also some scripts that are useful for producing a graphical output in R.
These are:
-format_for_graph.py
	Generates a file that can be easily graphed.
-rewrite_with_contig_of_interest.py
	Keeps only the contigs/species that the user wants to graph from the above output. This requires the user 
	to create a text file with the list of contigs/species they wish to view, written in the same way they are in the
	above output. The file should be saved as species_to_graph.txt
-graph_contig_of_interest.R
	Produces a pdf of the graph for each contig/species the user wanted to graph.


Samples of these can be found in the Sample_Results_Graphical folder.
