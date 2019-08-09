package statisticMerger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;


/**
 * Default implementation of a "Statistic Merger"
 * @author Michael Kluge
 *
 */
public abstract class AbstractStatisticMerger {
	
	public static final String FILE_ENDING = ".txt";
	protected static final LinkedHashMap<String, ArrayList<String>> MODULE_NAMES = new LinkedHashMap<String, ArrayList<String>>();
	
	/**
	 * Extracts the stuff from the single files and merges it into one file.
	 * @param inputPath
	 * @param outputPath
	 * @throws Exception
	 */
    protected void execute(String inputPath, String outputPath) throws Exception {
    	// create output folder
    	File outputFolder = new File(outputPath);
    	if(!outputFolder.isDirectory())
    		if(!outputFolder.mkdirs())
    			throw new IllegalArgumentException("Could not create output folder: '" + outputPath + "'!");
    	
    	// test if input folder is there
    	if(!new File(inputPath).isDirectory()) {
    		throw new IllegalArgumentException("Input path '" + inputPath + "' was not found.");
    	}
    		
    	// get all files in folder
    	ArrayList<String> dataFiles = this.getStatisticFiles(inputPath);
 
    	// run through all modules
    	for(String moduleName : this.getModuleNames()) {
	    	String moduleFileName = moduleName.replace(" ", "_") + FILE_ENDING;
	    	boolean writeHeader = true;
	    	String outfile = new File(outputFolder + File.separator + moduleFileName).getAbsolutePath();
	    	new File(outfile).delete();
	    	BufferedWriter outfileBW = new BufferedWriter(new FileWriter(outfile));
	    	
	    	// run though all files
	    	for(String file : dataFiles) {
	    		extractTable(moduleName, file, outfileBW, writeHeader);
	    		writeHeader = false;
	    	}
	    	
	    	outfileBW.flush();
	    	finalize(outfileBW);
	    	// close outfile
	    	outfileBW.close();
    	}
    }
	
	/**
	 * Gets all statistic files in a folder and in subfolders
	 * @param path Path to start with
	 * @return
	 */
	protected ArrayList<String> getStatisticFiles(String path) {
		return getStatisticFiles(path, new ArrayList<String>());
	}
	
	/**
	 * Gets all statistic files in a folder and in subfolders
	 * @param path Path to start with
	 * @param files Files which were already found
	 * @return
	 */
	private ArrayList<String> getStatisticFiles(String path, ArrayList<String> files) {
		File fpath = new File(path);

		// test all files in folder
		for(File f : fpath.listFiles()) {
			// if name is ok, add it
			if(f.isFile() && f.getName().matches(this.getNameOfStatisticFile())) {
				files.add(f.getAbsolutePath());
			}
			// call function recursive
			else if(f.isDirectory()) {
				getStatisticFiles(f.getAbsolutePath(), files);
			}
		}
		return files;
	}
	
	/**
	 * adds a new module names
	 * @param mergerName
	 * @param moduleNames
	 */
	public static void addModuleNames(String mergerName, ArrayList<String> moduleNames) {
		MODULE_NAMES.put(mergerName, moduleNames);
	}

	/**
	 * get the module names for a specific merger
	 * @param mergerName
	 * @return
	 */
    public static ArrayList<String> getModuleNames(String mergerName) {
    	return MODULE_NAMES.get(mergerName);
    }
    
    /**
     * Names of the modules
     * @return
     */
    public ArrayList<String> getModuleNames() {
    	return getModuleNames(this.getMergerName());
    }
	
    /****************************** ABSTRACT METHODS **********************************/
    
	/**
	 * Extracts the needed data out of the file for a specific module
	 * @param moduleName name of the module
	 * @param file Input file
	 * @param outfile output file
	 * @param writeHeader true, if header must be written
	 * @return
	 */
	protected abstract boolean extractTable(String moduleName, String file, BufferedWriter outfile, boolean writeHeader);
	
    
    /**
     * Name of the file where the statistic should be extracted
     * @return
     */
    public abstract String getNameOfStatisticFile();
    
    /**
     * Name of the statistic merger
     * @return
     */
    public abstract String getMergerName();
    
    /**
     * Is called after all files for a module were processed
     * @param outfile
     */
	public abstract void finalize(BufferedWriter outfile) throws IOException;
}