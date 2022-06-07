This folder contains the indices for aminoacids. Each subfolder contains the following:

- 000_raw_indices: The 57 indices as obtained from the ExPaSy ProtScale server on January 21st, 2019.
https://web.expasy.org/protscale/

- 001_processed_indices: The 57 indices after removing whitespace with the following command:

for infile in $(ls | rev | cut -c 5- | rev) 
do
	grep ':' ${infile}.txt | gsed -r 's/: /\t/g' > ../001_processed_indices/${infile}_processed.txt
done

- 002_all_data: The table with all the indices together. 

for infile in $(ls)
do
	diff <(cut -f 1 ${infile}) <(cut -f 1 coil_deleage_processed.txt)
done

I then used the following command to concatenate all the values:

paste <(cut -f 1 coil_deleage_processed.txt) <(cat ../002_all_data/aminoacid_one_letter_column.txt) > ../002_all_data/all_indices_table.txt

# Add the headers for the first two columns manually.
# With a loop, add each of the other indices to the table.

for infile in $(ls | rev | cut -c 5- | rev)
do
	cat <(echo ${infile}) <(cut -f 2 ${infile}.txt) > ../002_all_data/new_column.txt
	mv ../002_all_data/all_indices_table.txt ../002_all_data/tmp
	paste ../002_all_data/tmp ../002_all_data/new_column.txt > ../002_all_data/all_indices_table.txt
done

# Remove temporary files

# Remove extra spaces
sed 's/ //g' ../002_all_data/all_indices_table.txt > ../002_all_data/all_indices_final_table.txt





