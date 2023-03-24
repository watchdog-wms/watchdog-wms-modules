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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 *
 * @author friedel
 */
public class ATACSummarizer
{

    HashMap<String, ArrayList<PeakSummary>> summaries;
    ArrayList<Gene> genes;

    public ATACSummarizer()
    {
        this.summaries = new HashMap<>();
        this.genes = new ArrayList<>();
    }

    public void analyze(String infile)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));

            String[] content;

            String previousSample = null;
            String previousGene = null;

            String[] header = rd.readLine().split("\t");

            HashMap<String, Integer> nameToPos = new HashMap<>();
            ArrayList<PeakSummary> sampleSummaries = null;

            PeakSummary current = null;
            Gene gene;

            for (int i = 0; i < header.length; i++)
            {
                nameToPos.put(header[i], i);
            }

            String inline = rd.readLine();

            Peak peak;
            while (inline != null)
            {
                content = inline.split("\t");

                if (previousSample == null || !previousSample.equals(content[nameToPos.get("Sample")]))
                {

                    previousSample = content[nameToPos.get("Sample")];
                    sampleSummaries = new ArrayList<>();
                    this.summaries.put(previousSample, sampleSummaries);
                }
                if (previousGene == null || !previousGene.equals(content[nameToPos.get("geneId")]) || current == null)
                {
                    previousGene = content[nameToPos.get("geneId")];
                    current = null;
                    if (!previousGene.equals("NA"))
                    {
                        current = new PeakSummary();
                        gene = new Gene(previousGene, content[nameToPos.get("geneId")],
                                3 - 2 * Integer.parseInt(content[nameToPos.get("geneStrand")]),
                                content[nameToPos.get("seqnames")],
                                Integer.parseInt(content[nameToPos.get("geneStart")]), Integer.parseInt(content[nameToPos.get("geneEnd")]));

                        current.setGene(gene);

                        sampleSummaries.add(current);
                    }
                }
                if (current != null)
                {

                    int distanceToTSS = (int) Double.parseDouble(content[nameToPos.get("distanceToTSS")]);
                    /* if(content[nameToPos.get("distanceToTSS")].equals("1e+05"))
                    {
                        distanceToTSS=Integer.MAX_VALUE;
                        
                    }
                    else if(content[nameToPos.get("distanceToTSS")].equals("-1e+05"))
                    {
                        distanceToTSS=Integer.MIN_VALUE;
                        
                    }
                    else{
                        distanceToTSS=Integer.parseInt(content[nameToPos.get("distanceToTSS")]);
                    }*/

                    peak = new Peak(content[nameToPos.get("seqnames")],
                            Integer.parseInt(content[nameToPos.get("start")]), Integer.parseInt(content[nameToPos.get("end")]),
                            Double.parseDouble(content[nameToPos.get("V5")]),
                            content[nameToPos.get("V4")], content[nameToPos.get("annotation")], distanceToTSS);
                    current.addPeak(peak);
                }

                inline = rd.readLine();
            }
            rd.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ATACSummarizer.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void getUnassignedPeaks(String infile, String outfile, int minSize, int segmentMaxDist)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String previousSample = null;
            String previousGene = null;

            String[] header = rd.readLine().split("\t");

            HashMap<String, Integer> nameToPos = new HashMap<String, Integer>();
            ArrayList<PeakSummary> sampleSummaries = null;

            PeakSummary current = null;
            Gene gene;

            HashMap<String, Peak> peaks = new HashMap<String, Peak>();

            for (int i = 0; i < header.length; i++)
            {
                nameToPos.put(header[i], i);
            }

            String inline = rd.readLine();
            String[] content;
            Peak peak;
            HashMap<String, ArrayList<Peak>> peakPerAnno;
            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
            int countAll = 0;
            int countRemoved = 0;
            while (inline != null)
            {
                content = inline.split("\t");

                if (previousSample == null || !previousSample.equals(content[nameToPos.get("Sample")]))
                {

                    if (previousSample != null)
                    {
                        for (PeakSummary summary : this.summaries.get(previousSample))
                        {
                            peakPerAnno = summary.getPeaks();
                            for (String key : peakPerAnno.keySet())
                            {
                                for (Peak assignedPeak : peakPerAnno.get(key))
                                {

                                    if (peaks.containsKey(assignedPeak.getId() + "_" + previousSample))
                                    {
                                        peaks.remove(assignedPeak.getId() + "_" + previousSample);
                                        countRemoved++;
                                    }
                                }
                            }
                        }

                    }
                    System.out.println(previousSample + "\t" + countAll + "\t" + countRemoved + "\t" + (countAll - countRemoved));
                    previousSample = content[nameToPos.get("Sample")];

                    countAll = 0;
                    countRemoved = 0;
                }

                int distanceToTSS = (int) Double.parseDouble(content[nameToPos.get("distanceToTSS")]);
                peak = new Peak(content[nameToPos.get("seqnames")],
                        Integer.parseInt(content[nameToPos.get("start")]), Integer.parseInt(content[nameToPos.get("end")]),
                        Double.parseDouble(content[nameToPos.get("V5")]),
                        content[nameToPos.get("V4")] + "_" + previousSample, content[nameToPos.get("annotation")], distanceToTSS);

                if (peak.getAnnotation().equals("Distal Intergenic") && peak.getEnd() - peak.getStart() + 1 >= minSize)
                {
                    peaks.put(peak.getId(), peak);
                    countAll++;
                }

                inline = rd.readLine();
            }

            if (previousSample != null)
            {
                for (PeakSummary summary : this.summaries.get(previousSample))
                {
                    peakPerAnno = summary.getPeaks();
                    for (String key : peakPerAnno.keySet())
                    {
                        for (Peak assignedPeak : peakPerAnno.get(key))
                        {
                            if (peaks.containsKey(assignedPeak.getId() + "_" + previousSample))
                            {
                                peaks.remove(assignedPeak.getId() + "_" + previousSample);
                                countRemoved++;
                            }
                        }
                    }
                }
                System.out.println(previousSample + "\t" + countAll + "\t" + countRemoved + "\t" + (countAll - countRemoved));
            }

            Peak[] sorted_peaks = new Peak[peaks.size()];
            int i = 0;
            for (String key : peaks.keySet())
            {
                sorted_peaks[i] = peaks.get(key);
                i++;
            }
            Arrays.sort(sorted_peaks);

            int startSegment = 0;
            int maxEnd = sorted_peaks[startSegment].getEnd();
            System.out.println(sorted_peaks.length);
            wr.write("snapshotDirectory /mnt/raidproj/proj/projekte/expressionlab/RNATagging/Sequencing/HSV-1/Results_160217/analysis/IGV_Figures/\n");
            int countMock = 0;
            for (int j = 1; j < sorted_peaks.length; j++)
            {
                //  System.out.println(sorted_peaks[j].getId());
                if (!sorted_peaks[j].getChr().equals(sorted_peaks[j - 1].getChr()) || sorted_peaks[j].getStart() - maxEnd > segmentMaxDist)
                {

                    //  if(countMock==0)
                    {
                        wr.write("goto " + sorted_peaks[startSegment].getChr() + ":" + String.format(Locale.US, "%,d", sorted_peaks[startSegment].getStart() - segmentMaxDist) + "-" + String.format(Locale.US, "%,d", sorted_peaks[j - 1].getEnd() + segmentMaxDist) + "\n");
                        wr.write("snapshot " + sorted_peaks[startSegment].getChr() + "_" + sorted_peaks[startSegment].getStart() + "_" + sorted_peaks[j - 1].getEnd() + ".png\n");
                    }

                    startSegment = j;
                    maxEnd = sorted_peaks[startSegment].getEnd();
                    countMock = 0;
                } else
                {
                    if (sorted_peaks[j].getEnd() > maxEnd)
                    {
                        maxEnd = sorted_peaks[j].getEnd();
                    }
                }
                if (sorted_peaks[j].getId().contains("mock"))
                {
                    countMock++;
                }
            }
            if (startSegment == sorted_peaks.length - 1)
            {
                // if (countMock == 0)
                {
                    wr.write("goto " + sorted_peaks[startSegment].getChr() + ":" + String.format(Locale.US, "%,d", sorted_peaks[startSegment].getStart() - segmentMaxDist) + "-" + String.format(Locale.US, "%,d", sorted_peaks[startSegment].getEnd() + segmentMaxDist) + "\n");
                    wr.write("snapshot " + sorted_peaks[startSegment].getChr() + "_" + sorted_peaks[startSegment].getStart() + "_" + sorted_peaks[startSegment].getEnd() + ".png\n");
                }
            }

            wr.close();

            rd.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ATACSummarizer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void processPeaks(String infile, int dist1, int dist2)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String previousSample = null;
            String previousGene = null;

            String[] header = rd.readLine().split("\t");

            HashMap<String, Integer> nameToPos = new HashMap<>();
            ArrayList<PeakSummary> sampleSummaries = null;

            PeakSummary current = null;
            Gene gene;

            ArrayList<Peak> peaks = null;

            for (int i = 0; i < header.length; i++)
            {
                nameToPos.put(header[i], i);
            }

            String inline = rd.readLine();
            String[] content;
            Peak peak;
            while (inline != null)
            {
                content = inline.split("\t");

                if (previousSample == null || !previousSample.equals(content[nameToPos.get("Sample")]))
                {

                    if (previousSample != null)
                    {
                        System.out.println(previousSample);
                        sampleSummaries = assignPeaks(peaks, -1, dist1, dist2);
                        sampleSummaries.addAll(assignPeaks(peaks, 1, dist1, dist2));
                        this.summaries.put(previousSample, sampleSummaries);

                    }

                    previousSample = content[nameToPos.get("Sample")];
                    peaks = new ArrayList<Peak>();
                }

                int distanceToTSS = (int) Double.parseDouble(content[nameToPos.get("distanceToTSS")]);
                peak = new Peak(content[nameToPos.get("seqnames")],
                        Integer.parseInt(content[nameToPos.get("start")]), Integer.parseInt(content[nameToPos.get("end")]),
                        Double.parseDouble(content[nameToPos.get("V5")]),
                        content[nameToPos.get("V4")], content[nameToPos.get("annotation")], distanceToTSS);
                peaks.add(peak);

                inline = rd.readLine();
            }
            if (previousSample != null)
            {
                System.out.println(previousSample);
                sampleSummaries = assignPeaks(peaks, -1, dist1, dist2);
                sampleSummaries.addAll(assignPeaks(peaks, 1, dist1, dist2));
                this.summaries.put(previousSample, sampleSummaries);
            }

            rd.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ATACSummarizer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void processPeaksFromBED(String infile, String sample, int dist1, int dist2, boolean dOCRs)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));

            ArrayList<PeakSummary> sampleSummaries;

            ArrayList<Peak> peaks = null;

            String inline = rd.readLine();
            String[] content;
            Peak peak;
            String chr;
            int start;
            int end;
            double score;
            String id;
            peaks = new ArrayList<Peak>();
            while (inline != null)
            {
                content = inline.split("\t");

                chr = content[0];
                start = Integer.parseInt(content[1]);
                end = Integer.parseInt(content[2]);
                score = Double.parseDouble(content[4]);
                id = content[3];

                peak = new Peak(chr, start, end, score, id, null, -1);
                peaks.add(peak);

                inline = rd.readLine();
            }

            if (dOCRs)
            {
                sampleSummaries = assignPeaks(peaks, -1, dist1, dist2);
                sampleSummaries.addAll(assignPeaks(peaks, 1, dist1, dist2));
            } else
            {
              /*  sampleSummaries = assignPeaksInGenes(peaks, -1, dist1, dist2);
                sampleSummaries.addAll(assignPeaksInGenes(peaks, 1, dist1, dist2));*/
                
                sampleSummaries = assignPeaksInGenesAlternative(peaks, -1, dist1, dist2);
                sampleSummaries.addAll(assignPeaksInGenesAlternative(peaks, 1, dist1, dist2));
            }

            this.summaries.put(sample, sampleSummaries);

            rd.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ATACSummarizer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public ArrayList<PeakSummary> assignPeaksVersionPlosPathogens(ArrayList<Peak> peaks, int strand, int dist1, int dist2)
    {
        ArrayList<PeakSummary> result = new ArrayList<PeakSummary>();

        Peak[] sortedPeaks = peaks.toArray(new Peak[0]);
        Gene[] sortedGenes = this.genes.toArray(new Gene[0]);
        if (strand > 0)
        {
            Arrays.sort(sortedPeaks, new PeakSummary.ComparatorPositiveStrand());
            Arrays.sort(sortedGenes, new ComparatorPositiveStrand());
        } else
        {
            Arrays.sort(sortedPeaks, new PeakSummary.ComparatorNegativeStrand());
            Arrays.sort(sortedGenes, new ComparatorNegativeStrand());
        }

        int gene_pos = 0;
        int peak_pos = 0;
        int k = 0;
        int lastEnd = 0;

        PeakSummary peakSummary = null;

        while (gene_pos < sortedGenes.length)
        {
            if (sortedGenes[gene_pos].getStrand() == strand)
            {
                //System.out.println(sortedGenes[gene_pos].getSymbol()+"\t"+ sortedGenes[gene_pos].getChr()+"\t"+ sortedGenes[gene_pos].getEnd()+"\t"+ sortedPeaks[peak_pos]);

                peakSummary = new PeakSummary();
                peakSummary.setGene(sortedGenes[gene_pos]);

                //Find region
                while (peak_pos < sortedPeaks.length && sortedPeaks[peak_pos].getChr().compareTo(sortedGenes[gene_pos].getChr()) < 0)
                {
                    peak_pos++;
                }

                while (peak_pos < sortedPeaks.length && sortedPeaks[peak_pos].getChr().equals(sortedGenes[gene_pos].getChr())
                        && ((strand > 0 && sortedGenes[gene_pos].getEnd() > sortedPeaks[peak_pos].getEnd())
                        || (strand < 0 && sortedGenes[gene_pos].getStart() < sortedPeaks[peak_pos].getStart())))
                {
                    peak_pos++;

                }

                if (peak_pos < sortedPeaks.length)
                {
                    if ((strand > 0 && sortedGenes[gene_pos].getEnd() >= sortedPeaks[peak_pos].getStart() - dist1)
                            || strand < 0 && sortedGenes[gene_pos].getStart() <= sortedPeaks[peak_pos].getEnd() + dist1)
                    {

                        if (strand > 0)
                        {
                            lastEnd = sortedPeaks[peak_pos].getEnd();
                        } else
                        {
                            lastEnd = sortedPeaks[peak_pos].getStart();
                        }
                        k = 1;
                        peakSummary.addPeak(sortedPeaks[peak_pos]);

                        while (peak_pos + k < sortedPeaks.length && sortedPeaks[peak_pos + k].getChr().equals(sortedGenes[gene_pos].getChr())
                                && ((strand > 0 && ((sortedPeaks[peak_pos + k].getStart() - lastEnd + 1) <= dist2
                                || (sortedPeaks[peak_pos + k].getStart() - sortedGenes[gene_pos].getEnd() + 1) <= dist1))
                                || (strand < 0 && ((lastEnd - sortedPeaks[peak_pos + k].getEnd() + 1) <= dist2
                                || (sortedGenes[gene_pos].getStart() - sortedPeaks[peak_pos + k].getEnd() + 1) <= dist1)))) // da stand dist2, Fehler?
                        {
                            if (strand > 0)
                            {
                                lastEnd = sortedPeaks[peak_pos + k].getEnd();
                            } else
                            {
                                lastEnd = sortedPeaks[peak_pos + k].getStart();
                            }
                            peakSummary.addPeak(sortedPeaks[peak_pos + k]);
                            k++;
                        }
                    }
                }

                if (!peakSummary.isEmpty())
                {
                    result.add(peakSummary);
                }
            }
            gene_pos++;
        }

        return result;

    }

    public ArrayList<PeakSummary> assignPeaks(ArrayList<Peak> peaks, int strand, int dist1, int dist2)
    {
        ArrayList<PeakSummary> result = new ArrayList<>();

        Peak[] sortedPeaks = peaks.toArray(new Peak[0]);
        Gene[] sortedGenes = this.genes.toArray(new Gene[0]);
        if (strand > 0)
        {
            Arrays.sort(sortedPeaks, new PeakSummary.ComparatorPositiveStrand());
            Arrays.sort(sortedGenes, new ComparatorPositiveStrand());
        } else
        {
            Arrays.sort(sortedPeaks, new PeakSummary.ComparatorNegativeStrand());
            Arrays.sort(sortedGenes, new ComparatorNegativeStrand());
        }

        int gene_pos = 0;
        int peak_pos = 0;
        int k = 0;
        int lastEnd = 0;

        PeakSummary peakSummary = null;

        while (gene_pos < sortedGenes.length)
        {
            if (sortedGenes[gene_pos].getStrand() == strand)
            {
                //System.out.println(sortedGenes[gene_pos].getSymbol()+"\t"+ sortedGenes[gene_pos].getChr()+"\t"+ sortedGenes[gene_pos].getEnd()+"\t"+ sortedPeaks[peak_pos]);

                peakSummary = new PeakSummary();
                peakSummary.setGene(sortedGenes[gene_pos]);

                //Find region
                while (peak_pos < sortedPeaks.length && sortedPeaks[peak_pos].getChr().compareTo(sortedGenes[gene_pos].getChr()) < 0)
                {
                    peak_pos++;
                }

                while (peak_pos < sortedPeaks.length && sortedPeaks[peak_pos].getChr().equals(sortedGenes[gene_pos].getChr())
                        && ((strand > 0 && sortedGenes[gene_pos].getEnd() > sortedPeaks[peak_pos].getEnd())
                        || (strand < 0 && sortedGenes[gene_pos].getStart() < sortedPeaks[peak_pos].getStart())))
                {
                    peak_pos++;

                }

                if (peak_pos < sortedPeaks.length)
                {

                    if ((strand > 0 && sortedGenes[gene_pos].getEnd() >= sortedPeaks[peak_pos].getStart() - dist1)
                            || strand < 0 && sortedGenes[gene_pos].getStart() <= sortedPeaks[peak_pos].getEnd() + dist1)
                    {

                        if (strand > 0)
                        {
                            lastEnd = sortedPeaks[peak_pos].getEnd();
                        } else
                        {
                            lastEnd = sortedPeaks[peak_pos].getStart();
                        }
                        k = 1;
                        peakSummary.addPeak(sortedPeaks[peak_pos]);

                        while (peak_pos + k < sortedPeaks.length && sortedPeaks[peak_pos + k].getChr().equals(sortedGenes[gene_pos].getChr())
                                && ((strand > 0 && ((sortedPeaks[peak_pos + k].getStart() - lastEnd + 1) <= dist2
                                || (sortedPeaks[peak_pos + k].getStart() - sortedGenes[gene_pos].getEnd() + 1) <= dist1))
                                || (strand < 0 && ((lastEnd - sortedPeaks[peak_pos + k].getEnd() + 1) <= dist2
                                || (sortedGenes[gene_pos].getStart() - sortedPeaks[peak_pos + k].getEnd() + 1) <= dist1)))) // das stand dist2, Fehler?
                        {
                            if (strand > 0)
                            {
                                lastEnd = sortedPeaks[peak_pos + k].getEnd();
                            } else
                            {
                                lastEnd = sortedPeaks[peak_pos + k].getStart();
                            }
                            peakSummary.addPeak(sortedPeaks[peak_pos + k]);
                            k++;
                        }
                    }
                }

                if (!peakSummary.isEmpty())
                {
                    result.add(peakSummary);
                }
            }
            gene_pos++;
        }

        return result;

    }

    public ArrayList<PeakSummary> assignPeaksInGenes(ArrayList<Peak> peaks, int strand, int dist1, int dist2)
    {
        ArrayList<PeakSummary> result = new ArrayList<>();

        Peak[] sortedPeaks = peaks.toArray(new Peak[0]);
        Gene[] sortedGenes = this.genes.toArray(new Gene[0]);
        if (strand > 0)
        {
            Arrays.sort(sortedPeaks, new PeakSummary.ComparatorPositiveStrand());
            Arrays.sort(sortedGenes, new ComparatorPositiveStrandGeneStart());
        } else
        {
            Arrays.sort(sortedPeaks, new PeakSummary.ComparatorNegativeStrand());
            Arrays.sort(sortedGenes, new ComparatorNegativeStrandGeneStart());
        }

        int gene_pos = 0;
        int peak_pos = 0;
        int k = 0;
        int lastEnd = 0;

        PeakSummary peakSummary = null;

        while (gene_pos < sortedGenes.length)
        {

            if (sortedGenes[gene_pos].getStrand() == strand)
            {
                //System.out.println(sortedGenes[gene_pos].getSymbol()+"\t"+ sortedGenes[gene_pos].getChr()+"\t"+ sortedGenes[gene_pos].getEnd()+"\t"+ sortedPeaks[peak_pos]);

                peakSummary = new PeakSummary();
                peakSummary.setGene(sortedGenes[gene_pos]);

               

                //Find region
                while (peak_pos < sortedPeaks.length && sortedPeaks[peak_pos].getChr().compareTo(sortedGenes[gene_pos].getChr()) < 0)
                {
                    peak_pos++;

                }

                while (peak_pos < sortedPeaks.length && sortedPeaks[peak_pos].getChr().equals(sortedGenes[gene_pos].getChr())
                        && ((strand > 0 && sortedGenes[gene_pos].getStart() > sortedPeaks[peak_pos].getEnd())
                        || (strand < 0 && sortedGenes[gene_pos].getEnd() < sortedPeaks[peak_pos].getStart())))
                {
                    peak_pos++;

                }

                k = peak_pos;

                while (k < sortedPeaks.length && ((strand > 0 && sortedGenes[gene_pos].getEnd() > sortedPeaks[k].getStart())
                        || (strand < 0 && sortedGenes[gene_pos].getStart() < sortedPeaks[k].getEnd())))
                {

                    if ((sortedPeaks[k].getStart() >= sortedGenes[gene_pos].getStart() && sortedPeaks[k].getStart() <= sortedGenes[gene_pos].getEnd()) 
                     || (sortedPeaks[k].getEnd() >= sortedGenes[gene_pos].getStart() && sortedPeaks[k].getEnd() <= sortedGenes[gene_pos].getEnd()))
                    {
                        peakSummary.addPeak(sortedPeaks[k]);

                    }
                    k++;
                }

                if (!peakSummary.isEmpty())
                {
                    result.add(peakSummary);
                }
            }
            gene_pos++;
        }

        return result;

    }

    public ArrayList<PeakSummary> assignPeaksInGenesAlternative(ArrayList<Peak> peaks,  int strand, int dist1, int dist2)
    {
        ArrayList<PeakSummary> result = new ArrayList<>();

        Peak[] sortedPeaks = peaks.toArray(new Peak[0]);
        Gene[] sortedGenes = this.genes.toArray(new Gene[0]);
       
        Arrays.sort(sortedPeaks, new PeakSummary.ComparatorPositiveStrand());
        Arrays.sort(sortedGenes, new ComparatorPositiveStrandGeneStart());
        

        int gene_pos = 0;
        int peak_pos = 0;
        int k = 0;
        int lastEnd = 0;

        PeakSummary peakSummary = null;

        while (gene_pos < sortedGenes.length)
        {

            if (sortedGenes[gene_pos].getStrand() == strand)
            {
                

                peakSummary = new PeakSummary();
                peakSummary.setGene(sortedGenes[gene_pos]);

               

                //Find region
                while (peak_pos < sortedPeaks.length && sortedPeaks[peak_pos].getChr().compareTo(sortedGenes[gene_pos].getChr()) < 0)
                {
                    peak_pos++;

                }

                while (peak_pos < sortedPeaks.length && sortedPeaks[peak_pos].getChr().equals(sortedGenes[gene_pos].getChr()) &&
                         sortedPeaks[peak_pos].getEnd() < sortedGenes[gene_pos].getStart())
                {
                    peak_pos++;

                }

                k = peak_pos;

                while (k < sortedPeaks.length && sortedPeaks[k].getStart() < sortedGenes[gene_pos].getEnd())
                {

                    if ((sortedPeaks[k].getStart() <= sortedGenes[gene_pos].getStart() && sortedPeaks[k].getEnd()>= sortedGenes[gene_pos].getStart()) ||
                        (sortedPeaks[k].getStart() <= sortedGenes[gene_pos].getEnd() && sortedPeaks[k].getEnd()>= sortedGenes[gene_pos].getEnd()) 
                     || (sortedPeaks[k].getStart() >= sortedGenes[gene_pos].getStart() && sortedPeaks[k].getEnd()<= sortedGenes[gene_pos].getEnd()))
                    {
                        peakSummary.addPeak(sortedPeaks[k]);
                       
                    }
                    
                     
                    k++;
                }

                if (!peakSummary.isEmpty())
                {
                    result.add(peakSummary);
                }
            }
            gene_pos++;
        }

        return result;

    }
    
    public void readGenes(String infile)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));

            String[] content;

            String[] header = rd.readLine().split("\t");

            HashMap<String, Integer> nameToPos = new HashMap<>();
            ArrayList<PeakSummary> sampleSummaries = null;

            PeakSummary current = null;
            Gene gene;

            for (int i = 0; i < header.length; i++)
            {
                nameToPos.put(header[i], i);
            }

            String inline = rd.readLine();
            int strand;

            while (inline != null)
            {
                if (inline.startsWith("ENSG"))
                {
                    content = inline.split("\t");
                    strand = 1;
                    if (content[nameToPos.get("strand")].equals("-"))
                    {
                        strand = -1;
                    }

                    //wrong strand
                    //strand=strand*-1;
                    gene = new Gene(content[nameToPos.get("Geneid")], content[nameToPos.get("name")], strand,
                            content[nameToPos.get("chr")], Integer.parseInt(content[nameToPos.get("start")]),
                            Integer.parseInt(content[nameToPos.get("end")]));

                    this.genes.add(gene);
                }

                inline = rd.readLine();
            }
            rd.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ATACSummarizer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private class ComparatorPositiveStrand implements Comparator<Gene>
    {

        @Override
        public int compare(Gene o1, Gene o2)
        {
            if (o1.getChr().equals(o2.getChr()))
            {
                return o1.getEnd() - o2.getEnd();
            } else
            {
                return o1.getChr().compareTo(o2.getChr());
            }
        }

    }
    
    private class ComparatorPositiveStrandGeneStart implements Comparator<Gene>
    {

        @Override
        public int compare(Gene o1, Gene o2)
        {
            if (o1.getChr().equals(o2.getChr()))
            {
                return o1.getStart() - o2.getStart();
            } else
            {
                return o1.getChr().compareTo(o2.getChr());
            }
        }

    }

    private class ComparatorNegativeStrand implements Comparator<Gene>
    {

        @Override
        public int compare(Gene o1, Gene o2)
        {
            if (o1.getChr().equals(o2.getChr()))
            {
                return o2.getStart() - o1.getStart();
            } else
            {
                return o1.getChr().compareTo(o2.getChr());
            }
        }

        
    }
    
     private class ComparatorNegativeStrandGeneStart implements Comparator<Gene>
    {

        @Override
        public int compare(Gene o1, Gene o2)
        {
            if (o1.getChr().equals(o2.getChr()))
            {
                return o2.getEnd() - o1.getEnd();
            } else
            {
                return o1.getChr().compareTo(o2.getChr());
            }
        }

    }

    public void printToFile(String fileStart, int maxDistance, boolean countOverlapping) // für plos pathogens paper countOverlapping=false
    {
        try
        {

            for (String sample : this.summaries.keySet())
            {
                BufferedWriter wr = new BufferedWriter(new FileWriter(fileStart + "_" + sample + ".txt"));
                BufferedWriter wr2 = new BufferedWriter(new FileWriter(fileStart + "_" + sample + ".bed"));
                BufferedWriter wr3 = new BufferedWriter(new FileWriter(fileStart + "_" + sample + "_gene.bed"));
                wr.write("gene_id\tgene_symbol\tgene_chr\tgene_strand\tgene_start\tgene_end\tpromoter_id\t"
                        + "promoter_chr\tpromoter_start\tpromoter_end\tpromoter_score\tpromoter_width\tpromoter_annotation\tpromoter_distanceToTSS\t"
                        + "downstream_end\tdownstream_max\tdownstream_min\tdownstream_max_length\tdownstream_min_length\tdownstream_length\tdownstream_weightedlength\tdownstream_ids\n");

                ArrayList<PeakSummary> sampleSummaries = this.summaries.get(sample);

                for (PeakSummary summary : sampleSummaries)
                {
                    summary.summarize(maxDistance, countOverlapping);
                    wr.write(summary.toString() + "\n");
                    wr2.write(summary.getBEDDownstream());
                    wr3.write(summary.getDownStreamRegionAsBED());
                }

                wr.close();
                wr2.close();
                wr3.close();
            }
        } catch (IOException ex)
        {
            Logger.getLogger(ATACSummarizer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void printToFileInGene(String fileStart) // für plos pathogens paper countOverlapping=false
    {
        try
        {

            for (String sample : this.summaries.keySet())
            {
                BufferedWriter wr = new BufferedWriter(new FileWriter(fileStart + "_" + sample + ".txt"));
                BufferedWriter wr2 = new BufferedWriter(new FileWriter(fileStart + "_" + sample + ".bed"));
                wr.write("gene_id\tgene_symbol\tgene_chr\tgene_strand\tgene_start\tgene_end\tin_gene_length\tfraction_covered\tingene_ids\n");

                ArrayList<PeakSummary> sampleSummaries = this.summaries.get(sample);

                for (PeakSummary summary : sampleSummaries)
                {
                    summary.summarizeInGenes();
                    wr.write(summary.toStringInGenes() + "\n");
                    wr2.write(summary.getBEDInGene());
                }

                wr.close();
                wr2.close();
            }
        } catch (IOException ex)
        {
            Logger.getLogger(ATACSummarizer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void main(String[] args) throws ParseException
    {
        // create Options object
        Options options = new Options();

        options.addRequiredOption("b", "bed", true, "input BED file");
        options.addRequiredOption("o", "outputdir", true, "output directory");
        options.addRequiredOption("d1", "distance1", true, "maximum distance to gene end");
        options.addRequiredOption("d2", "distance2", true, "maximum distance to last added dOCR");
        options.addRequiredOption("n", "name", true, "Sample name");
        options.addRequiredOption("a", "annotation", true, "gene annotation");
        options.addOption("co", "countOverlapping", false, "count also dOCRs overlapping the gene");
        options.addOption("g", "gene", false, "count OCRs in genes");

        if (args.length == 0)
        {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("dOCR calculator", options);
            System.exit(0);
        }

        CommandLineParser parser = new DefaultParser();

        CommandLine cmd = parser.parse(options, args);

        boolean countOverlapping = false;

        String bedfile = cmd.getOptionValue("b");
        String outputdir = cmd.getOptionValue("o");
        int dist1 = Integer.parseInt(cmd.getOptionValue("d1"));
        int dist2 = Integer.parseInt(cmd.getOptionValue("d2"));
        String name = cmd.getOptionValue("n");
        String annotation = cmd.getOptionValue("a");

        if (cmd.hasOption("co"))
        {
            countOverlapping = true;
        }

        ATACSummarizer atac = new ATACSummarizer();
        atac.readGenes(annotation);

        atac.processPeaksFromBED(bedfile, name, dist1, dist2, true);

        atac.printToFile(outputdir + "Peak_assignment", 0, countOverlapping);

        if (cmd.hasOption("g"))
        {

            atac = new ATACSummarizer();
            atac.readGenes(annotation);
            atac.processPeaksFromBED(bedfile, name, dist1, dist2, false);
            atac.printToFileInGene(outputdir + "Gene_peak_assignment");
        }

    }
}
