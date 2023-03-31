
#first argument: exported expression matrix
#second argument: output file name (must be different than the first one)
cat utils/header.txt > $2
cat $1 >> $2
