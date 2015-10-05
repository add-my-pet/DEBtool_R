
#####	Function which "source" all the functions in a directory

## Argument
#
# path : the full name of the directory path

call_function <- function(path){

	# get the subdirectories path
	dir <- list.dirs(path, full.names=T)
	
	# Source all the functions in each directory
	for (i in 1:length(dir)){

		fin <- list.files(dir[i], full.names=T, include.dirs=F) # get functions path
		
		for (j in 1:length(fin)){
				source(fin[j])
		}	
	}
}


#######################
##### Test the function

path <- "/Users/gigamac/Desktop/FUN" # one directory with several subdirectories with n functions included

call_function(path=path)


