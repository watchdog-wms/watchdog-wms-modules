package statisticMerger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class CutadaptStatisticMerger extends AbstractStatisticMerger {
	
	public static final String TAB = "\t";
	public static final String HEADER = "parameter" + TAB + "value" + TAB + "position" + TAB + "sample";
	public static final String ENDING = "\\.fastq\\.cutadapt\\.out";
	public static final String DATA_FILE = ".*" + ENDING + "$";
	protected static final String MERGER_NAME = "CutadaptMerger";
	public static final String ADPATER_SEP = "=== Adapter [0-9]+ ===";
	public static final String SEQUENCE_START = "Sequence: ";
	public static final String PRIM5 = " times, it overlapped the 5' end of a read";
	public static final String PRIM3 = " times, it overlapped the 3' end or was within the read";
	public static final String COMMAND = "Command";
	
	// Available module names
	private static final String MODULE_CUTADAPT = "Cutadapt";
	protected static final ArrayList<String> MODULE_NAMES = new ArrayList<>();
	private static final ArrayList<String> ENDS = new ArrayList<>();
	
	 /**
     * add the used settings
     */
    static {
    	ENDS.add("5");
    	ENDS.add("3");
    	
        // add modules
        MODULE_NAMES.add(MODULE_CUTADAPT);
        
        addModuleNames(MERGER_NAME, MODULE_NAMES);
    }
	
	@Override
	protected boolean extractTable(String moduleName, String file, BufferedWriter outfile, boolean writeHeader) {
		File f = new File(file);
		// test if correct file is given
		if(!f.isFile() || !f.canRead() || !f.getName().matches(DATA_FILE))
			return false;
		
		// open file
		try {
			// write header
			if(writeHeader) {
				outfile.write(HEADER);
				outfile.newLine();
			}
			
			// process the file
			String line;
			BufferedReader r = new BufferedReader(new FileReader(f));

			boolean preAdapterSection = true;
			String adapter = "";
			String end = "";
			int removed = 0;
			
			String name = f.getName().replaceFirst(ENDING, "").replaceAll("^.", "");
			// read lines
			while((line = r.readLine()) != null) {
				line = new String(line.getBytes("US-ASCII"));

				if(preAdapterSection) {
					// check, if now adapter details are given.
					if(line.matches(ADPATER_SEP))
						preAdapterSection = false;
					else {
						String tmp[] = line.split(":");
						// check, if only one ":" was found
						if(tmp.length == 2) {
							String v = tmp[1].trim();
							
							if(!tmp[0].trim().startsWith(COMMAND))
								v = v.split(" ")[0];

							outfile.write(tmp[0].trim());
							outfile.write(TAB);
							outfile.write(v);
							outfile.write(TAB);
							outfile.write("0");
							outfile.write(TAB);
							outfile.write(name);
							outfile.newLine();
						}
					}
				}
				// get information about each of the adapters
				else { 
					if(line.startsWith(SEQUENCE_START)) {
						// get the adapter
						adapter = line.split(";")[0].replaceFirst(SEQUENCE_START, "").split(" ")[0];
						// get number of trimmed reads
						for(int i = 0; i < 2; i++) {
							line = r.readLine();
							
							if(line != null) {
								if(line.endsWith(PRIM5)) {
									removed = Integer.parseInt(line.replace(PRIM5, ""));
									end = "5";
								}
								else if(line.endsWith(PRIM3)) {
									removed = Integer.parseInt(line.replace(PRIM3, ""));
									end = "3";
								}
								// add dummy values if it does not match
								else {
									removed = 0;
									end = ENDS.get(i);
								}
							}
							// add dummy values if it does not match
							else {
								removed = 0;
								end = ENDS.get(i);
							}
							outfile.write(adapter);
							outfile.write(TAB);
							outfile.write(Integer.toString(removed));
							outfile.write(TAB);
							outfile.write(end);
							outfile.write(TAB);
							outfile.write(name);
							outfile.newLine();
						}
					}
				}
			}
			outfile.flush();
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
	public void finalize(BufferedWriter outfile) throws IOException {}
}