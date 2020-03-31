#!/bin/bash
## makeGFF.sh



FILE_PATH=`dirname $0`
number_of_files="3"
gff_file="None"
while getopts "i:j:f:c:n:g:l:ph" arg; do
  case $arg in
    i)
      file1_input=$OPTARG
      #echo $file1
      ;; 
    j)
      fasta_input=$OPTARG
      #echo $file1
      ;; 
    f)
      file_type=$OPTARG
      #echo $output_file
      ;;        
    c)
      contig=$OPTARG
      #echo $output_file
      ;;  
    n)
      number_of_files=$OPTARG
      #echo $output_file
      ;;
    g)
      gff_file=$OPTARG
	  ;;
    l)
      genome_lookup=$OPTARG
	  ;;
  
    h)
echo 'makeGFF.sh: Takes the output from predvirushost and make gff files for'
echo 'each contig.'
echo 'Version 1.0 2020'
echo '# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo 'usage: [options] <filename> <file type>'
echo ' '
echo 'Basic options:'
echo '  -i : <filename> proteins.txt output from predvirushost.sh'
echo '  -j : <filename> fasta file'
echo '  -f : <file type> the format of the fasta file (MGRAST, PROKKA, Refseq, Other)'
echo '  -c : <contig name> The name of the contig to make GFF for'
echo '  -n : <Number of contigs> if -c is not set, select the top n (ordered by number of proteins) contigs (Default is 3)'
echo ' '
exit      
      ;;
     
  esac
done

if [[ -z $file1_input ]]; then
echo 'Error: proteins.txt file needed. Specify with -i <input file>'
echo ' '
echo 'Use -h for more help.'
echo ' '
exit
fi

if [[ -z $fasta_input ]]; then
echo 'Error: Fasta File Needed. Specify with -j <input file>'
echo ' '
echo 'Use -h for more help.'
echo ' '
exit
fi

if [[ -z $genome_lookup ]]; then
echo 'Error: genome_lookup file needed. Specify with -l <input file>'
echo ' '
echo 'Use -h for more help.'
echo ' '
exit
fi

FILE_PATH2=`pwd`

echo "making tmp.faa using $fasta_input"
cat $fasta_input | sed 's/ /_/g' | tr '\n' ' ' | sed 's/>/\n/g' | sed 's/ /\t/' | sed 's/ //g' > tmp.faa

echo "$FILE_PATH2"
mv tmp.faa $FILE_PATH2



if [[ -z $contig ]]; then
echo 'running makeGFF.R'
${FILE_PATH}/makeGFF.R -p $FILE_PATH -q $FILE_PATH2 -f $file_type -i $file1_input -l $genome_lookup -n $number_of_files -g $gff_file
else
echo "running makeGFF.R using $contig"

${FILE_PATH}/makeGFF.R -p $FILE_PATH -q $FILE_PATH2 -f $file_type -i $file1_input -l $genome_lookup -c $contig -g $gff_file
fi




echo "Done."
