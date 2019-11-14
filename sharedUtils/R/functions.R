# writes a file to let Watchdog that script run to end and end
endNice <- function(fileName, params = NULL) {
	# write a file that we know, the script run to its end
	if(!is.null(fileName)) {
		file.create(fileName, showWarnings = FALSE)
		if(!is.null(params)) {
			write.table(t(data.frame(params)), file = fileName, quote = F, col.names = F, sep="\t")
		}
		write("?EOF!", file = fileName, append=TRUE, sep = "\n")
	}
	quit('no')
}

# test if an R package is installed on the system
testPackage <- function(pname) {
	return(is.element(pname, installed.packages()[,1]))
}

# test if an list of R packages is installed on the system and if not exits
endErrorMissingPackage <- function(pnames) {
	missing <- FALSE
	for(pname in pnames) {
		if(!testPackage(pname)) {
			missing <- TRUE
			print(paste('[ERROR] R library ', pname, ' is not installed.', sep = ''))
		}
	}
	if(missing) {
		quit('no', 12) # missing tools
	}
}

# loads a list of packages
loadPackages <- function(pnames) {
	for(pname in pnames) {
		print(paste('[INFO] loading R library ', pname, '...', sep = ''))
		library(pname, character.only = TRUE)
	}
}

# tests if a parameter is not empty
verifyInputNotEmpty <- function(inputOpt, name, exit) {
	input <- inputOpt[name][[1]]
	if(is.null(input) || input == '') {
		print(paste('[ERROR] Parameter --', name, ' can not be empty!', sep = ''))
		if(exit) { quit('no', 11) } # missing arguments
		return(FALSE)
	}
	if(!exit) { return(TRUE) }
}

# tests if a parameter is not empty
verifyFileExistence <- function(fname, exit) {
	# test if exists
	if(!file.exists(fname)) {
		print(paste('[ERROR] File ', fname, ' does not exist!', sep = ''))
		if(exit) { quit('no', 15) } # missing file
		return(FALSE)
	}
	# test if readable
	if(!file.access(fname, mode = 4)[1] == 0) {
		print(paste('[ERROR] File ', fname, ' is not readable!', sep = ''))
		if(exit) { quit('no', 15) } # missing file
		return(FALSE)
	}
	if(!exit) { return(TRUE) }
}

# retuns a output filename
getOutputFile <- function(parent, filename) {
	return(paste(parent, '/', filename, sep = ''))
}
