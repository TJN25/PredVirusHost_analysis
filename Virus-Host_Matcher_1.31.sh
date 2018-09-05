#!/bin/bash
#This folder should be the directory and have the sequence of interest in it in the form of name.faa
#This calls other scripts, makes some simple adjustments to the files and sorts the final results into a folder.
#See each script for details on what they do or look at the ReadMe.txt file.
#./Virus-Host_Matcher_1.31.sh -i test.faa -d _ -s 0 -e 3 -t y -r n -c 1


contig_cutoff=5
keep_files='n'
sig_hmms_only='y'



while getopts "i:d:s:e:t::r::c::" arg; do
  case $arg in
    i)
      file1=$OPTARG
      #echo $file1
      ;;
    d)
      delim=$OPTARG
      #echo $delim
      ;;
    s)
      pos_start=$OPTARG
      #echo $pos_start
      ;;
    e)
      pos_end=$OPTARG
      #echo $pos_end
      ;;      
    t)
      keep_files=$OPTARG
      #echo $keep_files
      ;;      
    r)
      sig_hmms_only=$OPTARG
      #echo $sig_hmms_only
      ;;      
    c)
      contig_cutoff=$OPTARG
      #echo $contig_cutoff
      ;;      
  esac
done




if [[ $contig_cutoff = '1' ]]; then
fasta_file=fastafile.txt
echo 'Using'
echo $fasta_file
echo 'as the input fasta file (fastafile.txt is the full file, fastafile_subset.txt is the smaller file)'
else
fasta_file=fastafile_subset.txt
echo 'Using'
echo $fasta_file
echo 'as the input fasta file (fastafile.txt is the full file, fastafile_subset.txt is the smaller file)'
fi


if [[ $sig_hmms_only = 'y' ]]; then
arVOG_hmm=arVOG_sig.hmm
baPOG_hmm=baPOG_sig.hmm
euVOG_hmm=euVOG_sig.hmm
echo 'Using Significant Models'
else
arVOG_hmm=arVOG-all.hmm
baPOG_hmm=baPOG_all.hmm
euVOG_hmm=euVOG_all.hmm
echo 'Using All Models'
fi

if [[ -z $file1 ]]; then
echo 'Error: Fasta File Needed'

elif [ $file1 = '-h' ]; then
echo 'Virus-Host_Matcher.sh: compares proteins from a number of contigs to hmm models and then scores'
echo 'each contig based on similarity to those models.'
echo 'Version 1.0 2016'
echo 'usage: <filename> <delimiter> <start_pos> <end_pos>'
echo '<filename> should be in a protein fasta format and include information about the contig or species'
echo '<delim> is the character that separates the contig name from the protein name/detail (if this is | then type __)'
echo '<start_pos> is the number of delimiters that you must pass to be left with just the contig name.'
echo '<end_pos> is the number of delimiters from the end needed to be left with just the contig'
echo 'For MG-RAST files the format would be <filename> _ 0 3 5 n y' 
echo 'which would keep all contigs with 5 or more proteins and remove the temp files once completed and would remove' 
echo 'some of the less significant models'

else
 cat $file1 > fastafile.txt #needed for running other scripts
chmod +x format_fastafile_to_remove_small_contigs.py
./format_fastafile_to_remove_small_contigs.py
 chmod +x remove_small_contigs.py 
 ./remove_small_contigs.py $delim $pos_start $pos_end $contig_cutoff
 echo 'Match '$file1' to HMM virus models'
 echo $arVOG_hmm
 echo 'Running HMMSearch on '$file1' and arVOG models' 
# each of the hmmsearch commands compare the input file '$file1' to one of the HMM model datasets
# the renaming is done to enable the files to be mre easily used with subsequent scripts
 hmmsearch --tblout arVOG_file.txt --noali $arVOG_hmm $fasta_file > $file1.output_ar.txt 
 echo 'Running HMMSearch on '$file1' and euVOG models' 
 hmmsearch --tblout euVOG_file.txt --noali  $euVOG_hmm $fasta_file > $file1.output_eu.txt
 echo 'Running HMMSearch on '$file1' and baPOG models' 
 hmmsearch --tblout baPOG_file.txt --noali $baPOG_hmm $fasta_file > $file1.output_ba.txt
 echo 'Making Protein List'
# this obtains a list of all the proteins in the input file
 grep '>' $fasta_file | cut -d ' ' -f1 | cut -d '>' -f2 > proteinnames_list.txt

chmod +x match_proteins_to_contigs.py
#this step match the protein to a genome/contig which allows for an overall score for each contig/genome to
#be calculated if more than one has been added.
./match_proteins_to_contigs.py $delim $pos_start $pos_end
echo 'Writing results into file'
chmod +x find_and_replace_space_with_tab.py
# the HMMsearch output includes a number of spaces and to work with the file more easily these are changed to one tab
./find_and_replace_space_with_tab.py
# only the protein name, e value and domain are kept from the output
grep 'arVOG_' arVOG_file.tab | cut -f1-3 > arVOG_list.txt
grep 'euVOG_' euVOG_file.tab | cut -f1-3 > euVOG_list.txt
grep 'bactPOG_' baPOG_file.tab | cut -f1-3 > baPOG_list.txt
chmod +x generate_output_file.py
# a list of the protein, the lowest e value and the domain this e value is from are listed.
./generate_output_file.py
sort output_tmp.txt | uniq > output.txt
chmod +x significant_matches_output.py
# a list that includes only those proteins that a match was found for
./significant_matches_output.py
chmod +x output_with_contig_names.py
# the complete list, including the species name for each protein is included
./output_with_contig_names.py


echo 'Calculating Scores'

chmod +x format_for_scoring.py
# This step assigns a minimun e value to the e values that are below the threshold in order to ensure
#that the score calculation is not weighhted to much based on the very small e values.
./format_for_scoring.py
# Calculates a score indicating which domain the genome likely belongs to
Rscript score_matches_for_contigs.R

#chmod +x format_output_scores.py
#./format_output_scores.py

keep_files_=$keep_files

fi

# these steps are needed to rename the files and put them in a folder

if [[ $keep_files_ = 'y' ]]; then

mkdir $file1.temp.folder
mv output.txt output_tmp.txt output_sig_hits_only.txt temp_file.txt scores_calc.txt output_scores.txt proteins_and_contigs_names.txt fastafile.txt output_with_species.txt $file1.temp.folder
mv arVOG_file.tab euVOG_file.tab baPOG_file.tab $file1.temp.folder
mv arVOG_file.txt baPOG_file.txt euVOG_file.txt $file1.temp.folder
mv arVOG_list.txt baPOG_list.txt euVOG_list.txt $file1.temp.folder
mv proteinnames_list.txt fastafile_subset.txt $file1.temp.folder
mv $file1.output_ar.txt $file1.output_ba.txt $file1.output_eu.txt $file1.temp.folder
command cat $file1.temp.folder/output.txt > $file1.txt
command cat $file1.temp.folder/output_sig_hits_only.txt > $file1.sig_hits_only.txt
command cat $file1.temp.folder/output_with_species.txt > $file1.output_with_species.txt
command cat $file1.temp.folder/output_scores.txt > $file1.scores.txt
command mkdir results.$file1.folder
command mv $file1.* results.$file1.folder
echo 'Temp files being moved to temp folder'
echo 'Finished'

elif [[ $keep_files_ = 'n' ]]; then


mkdir $file1.temp.folder
mv archaeal_contigs.txt output.txt output_tmp.txt output_sig_hits_only.txt temp_file.txt scores_calc.txt output_scores.txt proteins_and_contigs_names.txt fastafile.txt output_with_species.txt $file1.temp.folder
mv arVOG_file.tab euVOG_file.tab baPOG_file.tab $file1.temp.folder
mv arVOG_file.txt baPOG_file.txt euVOG_file.txt $file1.temp.folder
mv arVOG_list.txt baPOG_list.txt euVOG_list.txt $file1.temp.folder
mv proteinnames_list.txt fastafile_subset.txt $file1.temp.folder
mv $file1.output_ar.txt $file1.output_ba.txt $file1.output_eu.txt $file1.temp.folder
command cat $file1.temp.folder/output.txt > $file1.txt
command cat $file1.temp.folder/output_sig_hits_only.txt > $file1.sig_hits_only.txt
command cat $file1.temp.folder/output_with_species.txt > $file1.output_with_species.txt
command cat $file1.temp.folder/output_scores.txt > $file1.scores.txt
command cat $file1.temp.folder/archaeal_contigs.txt > $file1.archaeal_contigs.txt
command mkdir results.$file1.folder
command mv $file1.* results.$file1.folder
rm results.$file1.folder/$file1.temp.folder/*
rmdir results.$file1.folder/$file1.temp.folder/

echo 'Temp files being moved to temp folder'
echo 'Finished'

else

echo "Did not complete the commands: If you wish to debug then include 'y' in the command."

fi