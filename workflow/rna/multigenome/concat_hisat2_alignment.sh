## Save command line arguments
folder=$1
output=$2

## Create final output with header
echo -e "sample\tpercent_aligned" > "$output"

## Loop over files in the input folder
for log in $folder/*; do
    # Extract filename without path
    sample=$(basename "$log" .log)
    
    # Extract alignment percent
    percent=$(grep "overall alignment" "$log" | sed 's/% overall alignment rate//g')
    
    # Append to output file
    echo -e "${sample}\t${percent}" >> "$output"
done