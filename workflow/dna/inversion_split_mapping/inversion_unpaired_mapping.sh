############################################################
############################################################
#          Identifying Inversions from Alignment           #
############################################################
############################################################


# Creator: Zethus Woon Avery
# Function: Identify paired reads that map separately from a BAM file within inversion breakpoints
# Output: A folder with all of the reads that map to an inversion breakpoint and their pair maps elsewhere 
# as well as a summary file that counts the number of reads per file input


############################################################
#       Format Arguments and Instructional Messages        #
############################################################

## create help message
Help()
{
   # Display Help
   echo
   echo "Usage: Identify inversions from unpaired alignments in a given region."
   echo
   echo "Syntax: inversion_unpaired_mapping.sh [h|b|i or f|o|e|k|h]"
   echo
   echo "Required Input:"
   echo "-b PATH              Folder containing BAM files."
   echo
   echo "Options:"
   echo "-i CHROM:START-END   Inversion location [default: no location]"
   echo "-f PATH              Path to file listing inversion beakpoints, one per line [default: no location]"
   echo "-o PATH              Folder name for output [default: ./output/]"
   echo "-e STRING            Extension trailing the BAM files [default: .bam]"
   echo "-k TRUE/FALSE        Keep filtered BAM files [default: FALSE]"
   echo
   echo "Help Message:"
   echo "-h                   Print this help message."
   echo
}

## initialize arguments and set default output
folder=""
inversion=""
inversion_file=""
output="output"
extension=".bam"
keep="FALSE"

## load arguments
while getopts ":h:f:i:o:e:b:?:" option; do
   case $option in
      h) # display Help
         Help
         exit
         ;;
      b) # save folder path
         folder=$OPTARG
         ;;
      i) # save inversion location
         inversion=$OPTARG
         ;;
      f) # save inversion location
         inversion_file=$OPTARG
         ;;
      o) # save output path
         output=$OPTARG
         ;;
      e) # save extension name
         extension=$OPTARG
         ;;
      k) # specify if BAMs should be saved
         keep=$OPTARG
         ;;
      \?) ## display message if invalid argument used
         echo "Error: Invalid argument used"
         Help
         exit
         ;;
   esac
done

## create error message if essential arguments omitted
if [ -z "$folder" ]; then
  echo "Error: Path to folder containing input BAMs is required."
  Help
  exit 1
fi



############################################################
#          Identifying Inversions from Alignment           #
############################################################

## silence error messages
exec 2> /dev/null

## create output folder
mkdir -p $output/bams
echo "Initializing Environment..."

## create output summary file
echo "sample\tcount" > $output/split_mapping_sample_count.tsv

## remove bams in specified
if [ $keep == "FALSE" ]; then
   ## identify read pairs where one read maps to a given location and the other maps elsewhere
   if ["$inversion" == "" && "$inversion_file" == ""]; then # if no inversion location is given, run for the entire file
      for i in $folder/*.bam; do # iterate over files in folder
         sample=$(basename $i $extension) # get basename for individual file
         echo "Extracting split reads from $i...\n"
         samtools view -h -F 14 $i | awk '$0 ~ /^@/ || $7 != "*"' | samtools view -b - > $output/bams/$sample.bam # filter for requirement and keep sam header
         samtools view -c $output/bams/$sample > $output/count.txt # count the number of reads in bam
         awk -v sample="$sample" 'BEGIN { FS=OFS="\t" } {print $0, sample}' $output/count.txt > $output/temp.txt # create temporary text file with sample name and read count
         cat $output/temp.txt >> $output/split_mapping_sample_count.tsv # add sample information to summary file
         rm $output/count.txt $output/temp.txt $output/bams/$sample.bam # remove intermediate files
      done
   elif ["$inversion" == ""]; then # if inversion file with multiple breakpoints is given, run for each breakpoint
      echo "sample\tcount\tinversion" > $output/split_mapping_sample_count.tsv
      for j in $(cat $inversion_file); do # save each breakpoint to run through
         rm -rf $output/bams
         mkdir -p $output/$j # create folder for each breakpoint
         for i in $folder/*.bam; do
            sample=$(basename $i $extension)
            echo "Extracting split reads from $i for inversion $j..."
            samtools view -h -F 14 $i $j | awk '$0 ~ /^@/ || $7 != "*"' | samtools view -b - > $output/$j/$sample.bam
            samtools view -c $output/$j/$sample.bam > $output/count.txt
            awk -v sample="$sample" -v inversion="$j" 'BEGIN { FS=OFS="\t" } {print $0, sample, inversion}' $output/count.txt > $output/temp.txt 
            cat $output/temp.txt >> $output/split_mapping_sample_count.tsv
            rm $output/count.txt $output/temp.txt $output/$j/$sample.bam
         done
         rm -rf $output/$j/
      done
   else # specify single inversion location if given
      for i in $folder/*.bam; do
         sample=$(basename $i $extension)
         echo "Extracting split reads from $i..."
         samtools view -h -F 14 $i $inversion | awk '$0 ~ /^@/ || $7 != "*"' | samtools view -b - > $output/$sample
         samtools view -c $output/bams/$sample > $output/count.txt
         awk -v sample="$sample" 'BEGIN { FS=OFS="\t" } {print $0, sample}' $output/count.txt > $output/temp.txt
         cat $output/temp.txt >> $output/split_mapping_sample_count.tsv
         rm $output/count.txt $output/temp.txt $output/bams/$sample.bam
      done
   fi
elif [ $keep == "TRUE" ]; then # do not remove BAMs
   if ["$inversion" == "" && "$inversion_file" == ""]; then
      for i in $folder/*.bam; do
         sample=$(basename $i $extension)
         echo "Extracting split reads from $i..."
         samtools view -h -F 14 $i | awk '$0 ~ /^@/ || $7 != "*"' | samtools view -b - > $output/bams/$sample.bam 
         samtools view -c $output/bams/$sample > $output/count.txt
         awk -v sample="$sample" 'BEGIN { FS=OFS="\t" } {print $0, sample}' $output/count.txt > $output/temp.txt
         cat $output/temp.txt >> $output/split_mapping_sample_count.tsv
         rm $output/count.txt $output/temp.txt
      done
   elif ["$inversion" == ""]; then
      echo "sample\tcount\tinversion" > $output/split_mapping_sample_count.tsv
      for j in $(cat $inversion_file); do
         rm -rf $output/bams
         mkdir -p $output/bams/$j
         for i in $folder/*.bam; do
            sample=$(basename $i $extension)
            echo "Extracting split reads from $i for inversion $j..."
            samtools view -h -F 14 $i $j | awk '$0 ~ /^@/ || $7 != "*"' | samtools view -b - > $output/$j/$sample.bam
            samtools view -c $output/$j/$sample.bam > $output/count.txt
            awk -v sample="$sample" -v inversion="$j" 'BEGIN { FS=OFS="\t" } {print $0, sample, inversion}' $output/count.txt > $output/temp.txt
            cat $output/temp.txt >> $output/split_mapping_sample_count.tsv
            rm $output/count.txt $output/temp.txt
         done
      done
   else
      for i in $folder/*.bam; do
         sample=$(basename $i $extension)
         echo "Extracting split reads from $i...\n"
         samtools view -h -F 14 $i $inversion | awk '$0 ~ /^@/ || $7 != "*"' | samtools view -b - > $output/$sample
         samtools view -c $output/bams/$sample > $output/count.txt
         awk -v sample="$sample" 'BEGIN { FS=OFS="\t" } {print $0, sample}' $output/count.txt > $output/temp.txt
         cat $output/temp.txt >> $output/split_mapping_sample_count.tsv
         rm $output/count.txt $output/temp.txt
      done
   fi
else
  echo "Error: Improper usage of arguments."
  Help
  exit 1
fi