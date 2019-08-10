package statisticMerger;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;

import com.beust.jcommander.JCommander;

/**
 * Runner for the statistic mergers
 * @author Michael Kluge
 *
 */
public class StatisticMergerRunner {
	
	public static final String FASTQC_MERGER = "FastQC";
	public static final String STAR_MERGER = "Star";
	public static final String BAM_MERGER = "BamstatsMerger";
	public static final String CUTADAPT_MERGER = "CutadaptMerger";
	public static final String FEATURE_COUNTS_MERGER = "FeatureCounts";
	public static final String FLAGSTAT_MERGER = "FlagstatMerger";
	
	public static final int MAX_TIMEOUT = 30000; // wait a max of 30 seconds until file is created
	private static final String MERGED_TYPE = "mergedType";
	private static final String MERGED_FILE = "mergedFile";
	private static final String TAB = "\t";
	private static final String SEP = ":";
	private static final String EOF = "?EOF!";
	
    public static void main(String[] args) throws Exception {
    	Parameters params = new Parameters();
		JCommander parser = new JCommander(params, args);
		
		// display the help
		if(params.help) {
			parser.usage();
			System.exit(0);
		}
		else {
    		// get input parameters
    		String type = params.type;
    		String inputFolder = params.inputDir;
    		String outputFolder = params.outputDir;
    		
    		AbstractStatisticMerger m = null;
    		
    		// check, if type is valid and create class
    		if(FASTQC_MERGER.equals(type)) {
    			m = new FastqcStatisticMerger();
    		}
    		else if(STAR_MERGER.equals(type)) {
    			m = new StarStatisticMerger();
    		}
    		else if(BAM_MERGER.equals(type)) {
    			m = new BamstatsStatisticMerger();
    		}
    		else if(CUTADAPT_MERGER.equals(type)) {
    			m = new CutadaptStatisticMerger();
    		}
    		else if(FEATURE_COUNTS_MERGER.equals(type)) {
    			m = new FeatureCountsStatisticMerger();
    		}
    		else if(FLAGSTAT_MERGER.equals(type)) {
    			m = new FlagstatStatisticMerger();
    		}
    		else {
    			System.err.print("[ERROR] Not allowed type was set: '"+ type +"'!'");
    			System.exit(1);
    		}
    		
    		// execute the command
    		if(m != null) {
    			m.execute(inputFolder, outputFolder);
    			if(params.returnFilePath != null) {
    				File f = new File(params.returnFilePath);
    				if(!f.getParentFile().exists())
    					f.getParentFile().mkdirs();
    				
    				ArrayList<String> names = new ArrayList<>();
    				for(String moduleName : m.getModuleNames()) {
    					names.add(new File(outputFolder + File.separator + moduleName.replace(" ", "_") + AbstractStatisticMerger.FILE_ENDING).getAbsolutePath());
    				}
    				
    				FileWriter fw = new FileWriter(f);
    				fw.write(MERGED_TYPE + TAB + type);
    				fw.write(System.lineSeparator());
    				fw.write(MERGED_FILE + TAB + StringUtils.join(names, SEP));
    				fw.write(System.lineSeparator());
    				fw.write(EOF);
    				fw.flush();
    				fw.close();

    				// ensure that the file is really there
    				FileWatcher watcher = new FileWatcher(f.toPath(), MAX_TIMEOUT);
    				watcher.wasSuccessfull();
    			}
    		}
    		System.exit(0);
    	}
    }
}
