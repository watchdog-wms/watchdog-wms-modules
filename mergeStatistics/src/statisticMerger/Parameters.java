package statisticMerger;

import com.beust.jcommander.Parameter;

public class Parameters {
	
	@Parameter(names={"-type", "-t"}, description="type of the statistic merger (allowed values: FastQC, Star, BamstatsMerger, CutadaptMerger, FeatureCounts, FlagstatMerger)", required=true)
	protected String type;
	
	@Parameter(names={"-inputDir", "-i"}, description="path to input folder ", required=true)
	protected String inputDir;
	
	@Parameter(names={"-outputDir", "-o"}, description="path to output folder", required=true)
	protected String outputDir;
	
	@Parameter(names={"-returnFilePath", "-r"}, description="path to a file where the return parameters are written", required=false)
	protected String returnFilePath;
	
	@Parameter(names={"--help", "-h"}, description="print usage message and exit", help=true)
	protected boolean help = false;
}