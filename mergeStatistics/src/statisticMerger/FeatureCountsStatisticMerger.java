package statisticMerger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class FeatureCountsStatisticMerger extends AbstractStatisticMerger {
	
	public static final String END_FILENAME = ".counts.summary";
	public static final String DATA_FILE = ".*" + END_FILENAME;
	public static final String NEWLINE = System.lineSeparator();
	public static final String TAB = "\t";
	public static final String FILE_ENDING = ".txt";
	protected static final String NAME_OF_FILE = "Filename";
	protected static final String MERGER_NAME = "FeatureCounts";
	
	private final StringBuffer CONTENT = new StringBuffer();
	
	// Available module names
	private static final String MODULE_ASSIGNMENT = "Assignment statistic";
	protected static final ArrayList<String> MODULE_NAMES = new ArrayList<String>();
	
	
	 /**
     * add the used settings
     */
	static {
		// add module name
		ArrayList<String> modules = new ArrayList<String>();
		modules.add(MODULE_ASSIGNMENT);
		FastqcStatisticMerger.addModuleNames(MERGER_NAME, modules);
	}
	
	@Override
	protected boolean extractTable(String moduleName, String file, BufferedWriter outfile, boolean writeHeader) {
		File f = new File(file);
		// test if correct file is given
		if(!f.isFile() || !f.canRead() || !f.getName().matches(DATA_FILE))
			return false;
		
		// open file
		try {
			BufferedReader r = new BufferedReader(new FileReader(f));
			
			StringBuffer headerB = new StringBuffer();
			StringBuffer valuesB = new StringBuffer();
			
			String line = null;
			String name = null;
			String value = null;
			
			// read lines
			while((line = r.readLine()) != null) {
				// try to split line
				String tmp[] = line.split(TAB);
				if(tmp.length == 2) {
					// replace white spaces at start and end
					name = tmp[0].trim();
					value = tmp[1].trim();
										
					// if header should be written
					if(writeHeader) {
						headerB.append(name);
						headerB.append(TAB);
					}
					valuesB.append(value);
					valuesB.append(TAB);
				}
			}
			// add the filename
			headerB.append(NAME_OF_FILE);
			valuesB.append(f.getName().replaceFirst(END_FILENAME, ""));	

			// write to buffer
			if(writeHeader) 
				CONTENT.append(headerB.toString());
			
			CONTENT.append(NEWLINE);
			CONTENT.append(valuesB.toString());
			
			// close the file
			r.close();
			return true;
		}
		catch(Exception e) {
			e.printStackTrace(); 
		}
		return false;
	}

	@Override
	public String getNameOfStatisticFile() {
		return DATA_FILE;
	}

	@Override
	public String getMergerName() {
		return MERGER_NAME;
	}

	@Override
	public void finalize(BufferedWriter outfile) throws IOException {
		outfile.write(CONTENT.toString());
		outfile.flush();
		this.CONTENT.delete(0, this.CONTENT.length());
	}
}
