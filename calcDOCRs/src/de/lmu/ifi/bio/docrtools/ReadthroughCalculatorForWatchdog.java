/*
 * The MIT License
 *
 * Copyright 2023 friedel.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package de.lmu.ifi.bio.docrtools;

import java.util.HashSet;
import org.apache.commons.cli.*;


/**
 *
 * @author friedel
 */
public class ReadthroughCalculatorForWatchdog
{

    public static final String DEFAULTWINDOWLENGTH = "5000";
    public static final String DEFAULTOVERLAP = "25";
    public static final String DEFAULTSTRAND = "0";

    public static final String DEFAULTNORM = "2E9";

    public static void main(String[] args)
    {

        Options options = new Options();

        Option annoFile = new Option("a", "annotation", true, "annotation file path");
        annoFile.setRequired(true);
        options.addOption(annoFile);

        Option output = new Option("o", "output", true, "output file");
        output.setRequired(true);
        options.addOption(output);

        Option input = new Option("i", "input", true, "input file");
        input.setRequired(true);
        options.addOption(input);

        Option rcounts = new Option("g", "genecounts", true, "gene read count file");
        rcounts.setRequired(true);
        options.addOption(rcounts);

        Option strandspecificOption = new Option("s", "strandedness", true, "strandedness: 0=not strandspecific, 1=first read indicates strand, 2=second read indicates strand");
        strandspecificOption.setRequired(false);
        options.addOption(strandspecificOption);

        Option connectionFiles = new Option("c", "connections", true, "connection file paths");
        connectionFiles.setRequired(false);
        options.addOption(connectionFiles);

        Option readthrough = new Option("t", "readthroughLength", true, "length of downstream window in which read-through is calculated");
        readthrough.setRequired(false);
        options.addOption(readthrough);

        Option readin = new Option("n", "readinLength", true, "length of upstream window in which read-in is calculated");
        readin.setRequired(false);
        options.addOption(readin);

        Option overlap = new Option("v", "overlap", true, "minimum overlap for read-through/in window");
        overlap.setRequired(false);
        options.addOption(overlap);

        Option idxstats = new Option("x", "idxstats", true, "idxstats file with numbers of mapped reads per chromosome");
        idxstats.setRequired(false);
        options.addOption(idxstats);

        Option norm = new Option("m", "normFactor", true, "normalizing factor for FPKM calculation");
        norm.setRequired(false);
        options.addOption(norm);

        Option excludeChr = new Option("e", "exclude", true, "chromosomes to exclude from calculating total mapped reads, separated by ,");
        excludeChr.setRequired(false);
        options.addOption(excludeChr);

        Option dOCRFile = new Option("d", "dOCRFile", true, "file containing dOCR lengths");
        dOCRFile.setRequired(false);
        options.addOption(dOCRFile);

        Option windowLengthOption = new Option("w", "windowLength", true, "size of windows for dOCR transcription analysis");
        windowLengthOption.setRequired(false);
        options.addOption(windowLengthOption);
        
        Option excludeTypeOption = new Option("z", "excludeType", true, "excludeTypes");
        excludeTypeOption.setRequired(false);
        options.addOption(excludeTypeOption);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try
        {
            cmd = parser.parse(options, args);

            String annotation = cmd.getOptionValue("annotation");

            String outputFile = cmd.getOptionValue("output");

            String file = cmd.getOptionValue("input");

            String genecounts = cmd.getOptionValue("genecounts");

            ReadOutFromBamCalculator calc = new ReadOutFromBamCalculator();
            calc.readGenes(annotation);

            String connections = cmd.getOptionValue("connections");

            String strandedness = cmd.getOptionValue("strandedness", DEFAULTSTRAND);

            String idxFile = cmd.getOptionValue("idxstats");

            double normalizingFactors = Double.parseDouble(cmd.getOptionValue("normFactor", DEFAULTNORM));

            HashSet<String> excludeNames = new HashSet<>();
            String[] excluded = cmd.getOptionValue("exclude").split(",");
            for (String e : excluded)
            {
                excludeNames.add(e);
            }
            
            String excludedTypesString = cmd.getOptionValue("excludeType");
            HashSet<String> excludedTypes = new HashSet<>();
            if (excludedTypesString != null)
            {
                String[] excludedType = excludedTypesString.split(",");
                for (String e : excludedType)
                {
                    excludedTypes.add(e);
                }

            } else
            {
                
                    excludedTypes.add("pseudogene");
                    excludedTypes.add("snRNA");
            }
            
        

            String dOCRFilename = cmd.getOptionValue("dOCRFile");
            

            if (idxFile != null)
            {
                calc.getMillionMappedReads(idxFile, normalizingFactors, excludeNames);
            }

            if (dOCRFilename != null)
            {
                
                if(idxFile == null)
                {
                     System.err.println("idxstats file has to be provided");
                     System.exit(1);
                }

                try
                {
                    int windowLength=Integer.parseInt(cmd.getOptionValue("windowLength"));
                    calc.matchReadthroughToDOCRs(file, dOCRFilename, outputFile, windowLength, Integer.parseInt(cmd.getOptionValue("overlap", DEFAULTOVERLAP)), strandedness.equals("2"), !strandedness.equals(DEFAULTSTRAND));
                    
                } catch (NumberFormatException e)
                {
                    System.err.println(e.getMessage());

                }
            } 
            else
            {

                if (connections == null)
                {

                    String tmpneg = outputFile + ".negstrand.gff3";
                    String tmpnegdist = outputFile + ".negstrand.distances.gff3";
                    String tmpnegpairs = outputFile + ".negstrand.pairs.gff3";

                    String tmppos = outputFile + ".posstrand.gff3";
                    String tmpposdist = outputFile + ".posstrand.distances.gff3";
                    String tmppospairs = outputFile + ".posstrand.pairs.gff3";

                    UTRCalculator utrcalc = new UTRCalculator();
                    utrcalc.readGeneInformationFromGTF(annotation, -1, excludedTypes, "([0-9]+|X|Y|MT)");
                    utrcalc.getConnections(tmpneg, tmpnegdist, tmpnegpairs, 5000, 300, -1);

                    utrcalc = new UTRCalculator();
                    utrcalc.readGeneInformationFromGTF(annotation, 1, excludedTypes, "([0-9]+|X|Y|MT)");
                    utrcalc.getConnections(tmppos, tmpposdist, tmppospairs, 5000, 300, 1);

                    connections = tmpneg + "," + tmppos;

                }

                calc.readConnections(connections.split(","));

                calc.getGeneReadCounts(genecounts);
                try
                {
                    calc.getReadoutReadin(file, outputFile,
                            Integer.parseInt(cmd.getOptionValue("readinLength", DEFAULTWINDOWLENGTH)),
                            Integer.parseInt(cmd.getOptionValue("readthroughLength", DEFAULTWINDOWLENGTH)),
                            Integer.parseInt(cmd.getOptionValue("overlap", DEFAULTOVERLAP)), strandedness.equals("2"), !strandedness.equals(DEFAULTSTRAND));

                } catch (NumberFormatException e)
                {
                    System.err.println(e.getMessage());

                }
            }

        } catch (ParseException e)
        {
            System.out.println(e.getMessage());
            formatter.printHelp("read-through calculator", options);

        }

    }

}
