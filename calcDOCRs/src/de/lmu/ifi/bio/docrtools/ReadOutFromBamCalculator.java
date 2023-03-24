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

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 *
 * @author friedel
 */
public class ReadOutFromBamCalculator
{

    private HashMap<String, Gene> genes;
    private HashMap<String, Integer> distances;
    private HashMap<String, Double> millionMappedReads;

    public ReadOutFromBamCalculator()
    {
        this.genes = new HashMap<>();
        this.millionMappedReads= new HashMap<>();

    }
    
    public void readDistances(String infile)
    {
       
        try
        {
            distances= new HashMap<>();
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String inline = rd.readLine();
            String[] content;
            String gene;
            int distance;
            while (inline != null)
            {
                content = inline.split("\t");
                gene=content[0];
                distance=Integer.parseInt(content[8]);
                distances.put(gene, distance);
                
                inline = rd.readLine();
            }
            rd.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
       
    }

    public void readGenes(String gtfFile)
    {

        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(gtfFile));
            String inline = rd.readLine();
            String[] content;
            Gene gene;
            String biotype;

            String chr;
            int geneStart;
            int geneEnd;
            int geneStrand;
            String geneId;
            String geneName;

            String[] info;
            while (inline != null)
            {
                if (!inline.startsWith("#"))
                {
                    content = inline.split("\t");
                    if (content[2].equals("gene"))
                    {
                        chr = content[0];

                        if (chr.equals("chrMT"))
                        {
                            chr = "chrM";
                        }

                        geneStrand = 1;
                        if (content[6].equals("-"))
                        {
                            geneStrand = -1;
                        }
                        geneStart = Integer.parseInt(content[3]);
                        geneEnd = Integer.parseInt(content[4]);

                        info = content[8].split(";");

                        geneId = null;
                        geneName = null;
                        biotype = null;

                        for (String in : info)
                        {

                            if (in.startsWith("gene_id"))
                            {
                                geneId = in.replace("gene_id ", "").replace("\"", "").trim();
                            }
                            if (in.startsWith(" gene_name"))
                            {
                                geneName = in.replace("gene_name ", "").replace("\"", "").trim();
                            }
                            if (in.startsWith(" gene_biotype"))
                            {
                                biotype = in.replace("gene_biotype ", "").replace("\"", "").trim();
                            }
                        }

                        gene = new Gene(geneId, geneName, geneStrand, chr, geneStart, geneEnd);

                        if (!this.genes.containsKey(gene.getId()))
                        {
                            this.genes.put(gene.getId(), gene);
                        }
                    }

                }
                inline = rd.readLine();
            }
            rd.close();
        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
    
    public void getMillionMappedReads(String infile, double divisonFactor, HashSet<String> excludeNames)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String inline = rd.readLine();
            String[] content;
            
            String sample;
            String chr;
            double count;
 
            while (inline != null)
            {
                
                if(!inline.contains("\tmapped\t"))
                {
                    content = inline.split("\t");
                    chr=content[0];
                    count=Double.parseDouble(content[2]);
                    sample=content[4];
                    
                
                    
                    if(!chr.equals("*") && ( excludeNames== null || !excludeNames.contains(chr)))
                    {
                        if(!this.millionMappedReads.containsKey(sample))
                        {
                            this.millionMappedReads.put(sample, count);
                        }
                        else{
                            this.millionMappedReads.put(sample, count+this.millionMappedReads.get(sample));
                        }
                    }
                    
                }
                
                inline = rd.readLine();
            }
            rd.close();
            
            
            for(String key:this.millionMappedReads.keySet())
            {
                
                this.millionMappedReads.put(key, this.millionMappedReads.get(key)/divisonFactor);
               
            }
            
        } catch (IOException ex)
        {
            Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void getGeneReadCounts(String infile)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String inline = rd.readLine();
            String[] content;

            Gene gene;
            while (inline != null)
            {
                if (inline.startsWith("ENSG"))
                {
                    content = inline.split("\t");
                    gene = this.genes.get(content[0]);
                    if (gene != null)
                    {
                        gene.setExonLength(Integer.parseInt(content[5]));
                        gene.setExonReadCount(Integer.parseInt(content[6]));
                    }

                }

                inline = rd.readLine();
            }
            rd.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void matchReadthroughToDOCRs(String bamFile, String dOCRFile, String outfile,int windowlength, int minOverlap,boolean strandReversed, boolean strandSpecific)
    {
        
        try
        {
            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
            SamReaderFactory srf = SamReaderFactory.make();
            srf.validationStringency(ValidationStringency.SILENT);
            srf.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
            SamReader reader = srf.open(new File(bamFile));
            
            String name=bamFile.replace(".bam", "");
            name=name.replaceAll(".*/", "");

            BufferedReader rd = new BufferedReader(new FileReader(dOCRFile));
            
            
            String[] header=rd.readLine().split("\t");
            
            int indexEnd=-1;
            int indexLength=-1;
            for(int i=0;i<header.length;i++)
            {
                if(header[i].equals("downstream_end"))
                {
                    indexEnd=i;
                }
                if(header[i].equals("downstream_length"))
                {
                    indexLength=i;
                }
                
            }
            
            String chr;
            int strand;
            int start;
            int end;
            
            String inline = rd.readLine();
            String[] content;
            int windowStart;
            int windowEnd;
            int readCount;
            double normFactor=-1;
            
            if(this.millionMappedReads.containsKey(name))
            {
                   normFactor=this.millionMappedReads.get(name);

            }
            else{
                 System.err.println("sample "+name+" not included in idxstats file");
                 System.exit(1);
            }
            
            while (inline != null)
            {
              
                
                
                content = inline.split("\t");
             
                for(int i=0;i<6;i++)
                {
                    wr.write(content[i]+"\t");
                }
                 wr.write(content[indexEnd]+"\t"+content[indexLength]);
                
                chr=content[2];
                strand=Integer.parseInt(content[3]);
                
                
                if(strand>0)
                {
                    start=Integer.parseInt(content[5])+1;
                    end=Integer.parseInt(content[indexEnd]);
                    
                }
                else{
                    end=Integer.parseInt(content[4])-1;
                    start=Integer.parseInt(content[indexEnd]);
                    
                }
                
                
                windowStart=start;
                windowEnd=start+windowlength-1;
                
                if(strand<0)
                {
                    windowStart=end-windowlength+1;
                    windowEnd=end;
                }
                while((strand>0 & windowStart<end) || (strand<0 & windowEnd>start))
                {

                    
                    windowStart=windowStart+strand*windowlength;
                    windowEnd=windowEnd+strand*windowlength;
                    readCount = this.countReads(reader, chr, windowStart, windowEnd, minOverlap, strand, strandReversed, strandSpecific, false);
                   
                    double fpkm=(double)readCount/windowlength/normFactor;
                    wr.write("\t"+fpkm);
                
                    
                }

                wr.write("\n");
                inline = rd.readLine();
            }
            
            
            rd.close();
          
            wr.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
        } 
        
    }
    
    public void getReadoutReadin(String bamFile, String outfile, int readinLength, int readoutLength, int minOverlap, boolean strandReversed, boolean strandSpecific)
    {
        BufferedWriter wr = null;
        try
        {
            SamReaderFactory srf = SamReaderFactory.make();
            srf.validationStringency(ValidationStringency.SILENT);
            srf.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
            SamReader reader = srf.open(new File(bamFile));
            
            Gene gene;
            int startRegion;
            int endRegion;
            Gene connectedGene;
            int countReadIn;
            int countReadOut;
            wr = new BufferedWriter(new FileWriter(outfile));
            wr.write("Ensembl_ID\tChr\tStrand\tGene_expression\tReadin_start\tReadin_end\tReadin_count\tReadin\tReadout_start\tReadout_end\tReadout_count\tReadout\tDownstream_FPKM\n");

            double geneExpression;
            double readin;
            double readout;
            boolean verbose=false;
            
            String name=bamFile.replace(".bam", "");
            name=name.replaceAll(".*/", "");

            System.out.println(name);
                
            for (String g : genes.keySet())
            {
                gene = this.genes.get(g);

               
                connectedGene = gene.getUpstreamGene();

                countReadIn = -1;
                geneExpression = (double) gene.getExonReadCount() / gene.getExonLength();

                wr.write(gene.getId() + "\t" + gene.getChr() + "\t" + gene.getStrand() + "\t" + geneExpression);

                if (gene.getStrand() > 0)
                {
                    endRegion = gene.getStart() - 1;

                    startRegion = gene.getStart() - readinLength;

                    if (connectedGene != null)
                    {
                        if (startRegion > 0)
                        {
                            if (connectedGene.getEnd() >= startRegion)
                            {
                                startRegion = -1;
                            }
                        }
                    }

                    if (startRegion > 0)
                    {

                        countReadIn = this.countReads(reader, gene.getChr(), startRegion, endRegion, minOverlap, gene.getStrand(), strandReversed, strandSpecific, verbose);

                    }

                } else
                {

                    startRegion = gene.getEnd() + 1;

                    endRegion = gene.getEnd() + readinLength;

                    if (connectedGene != null)
                    {

                        if (connectedGene.getStart() <= endRegion)
                        {
                            startRegion = -1;
                        }

                    }

                    if (startRegion > 0)
                    {

                        countReadIn = this.countReads(reader, gene.getChr(), startRegion, endRegion, minOverlap, gene.getStrand(), strandReversed, strandSpecific, verbose);

                    }

                }

                wr.write("\t" + startRegion + "\t" + endRegion + "\t" + countReadIn);

                readin = (double) countReadIn / readinLength;
                readin /= geneExpression;

                if (countReadIn < 0 || geneExpression == 0 || (gene.getUpstreamGene() == null && !gene.isFirstInChromosome()))
                {
                    wr.write("\tNA");
                } else
                {
                    wr.write("\t" + readin);
                }

                
                connectedGene = gene.getDownstreamGene();

                countReadOut = -1;
                if (gene.getStrand() < 0)
                {
                    endRegion = gene.getStart() - 1;

                    startRegion = gene.getStart() - readoutLength;

                    if (connectedGene != null)
                    {
                        if (startRegion > 0)
                        {
                            if (connectedGene.getEnd() >= startRegion)
                            {
                                startRegion = -1;
                            }
                        }
                    }

                    if (startRegion > 0)
                    {

                        countReadOut = this.countReads(reader, gene.getChr(), startRegion, endRegion, minOverlap, gene.getStrand(), strandReversed, strandSpecific, verbose);

                    }

                } else
                {

                    startRegion = gene.getEnd() + 1;

                    endRegion = gene.getEnd() + readoutLength;

                    if (connectedGene != null)
                    {

                        if (connectedGene.getStart() <= endRegion)
                        {
                            startRegion = -1;
                        }

                    }

                    if (startRegion > 0)
                    {

                        countReadOut = this.countReads(reader, gene.getChr(), startRegion, endRegion, minOverlap, gene.getStrand(), strandReversed, strandSpecific, verbose);

                    }

                }
                
                
                wr.write("\t" + startRegion + "\t" + endRegion + "\t" + countReadOut);

                readout = (double) countReadOut / readoutLength;
                readout /= geneExpression;
                if (countReadOut < 0 || geneExpression == 0 || (gene.getDownstreamGene() == null && !gene.isLastInChromosome()))
                {
                    wr.write("\tNA");
                } else
                {
                    wr.write("\t" + readout);
                }
                
                
               
                
                if(this.millionMappedReads.containsKey(name) && countReadOut>=0)
                {
                    double downstream_FPKM=(double)countReadOut/readoutLength/this.millionMappedReads.get(name);
                    wr.write("\t"+downstream_FPKM);
                }
                else{
                    wr.write("\tNA");
                }
                
                wr.write("\n");

            }
            wr.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
        } finally
        {
            try
            {
                wr.close();
            } catch (IOException ex)
            {
                Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    public void getReadoutEnds(String bamFile, String outfile, int maxReadoutConsidered, boolean strandReversed, double threshold, int windowSize, boolean useGeneExpression)
    {
        BufferedWriter wr = null;
        BufferedWriter wrBed = null;
        System.out.println(bamFile);
        System.out.println(outfile);
        try
        {
            SamReaderFactory srf = SamReaderFactory.make();
            srf.validationStringency(ValidationStringency.SILENT);
            srf.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
            SamReader reader = srf.open(new File(bamFile));
            // get reads overlapping this region
            Gene gene;
            int startRegion;
            int endRegion;
            Gene connectedGene;
            int readoutEnd;
            wr = new BufferedWriter(new FileWriter(outfile));
            wrBed = new BufferedWriter(new FileWriter(outfile.replace(".txt", ".bed")));
            wr.write("Ensembl_ID\tChr\tStrand\tGene_expression\tRegion_start\tRegion_end\tReadout_end\n");

            double geneExpression;
            double readout;

            for (String g : genes.keySet())
            {
                gene = this.genes.get(g);

                // quantify readin;
                geneExpression = (double) gene.getExonReadCount() / gene.getExonLength();

                // quantify readout;
                connectedGene = gene.getDownstreamGene();

                readoutEnd = -1;
                if (gene.getStrand() < 0)
                {
                    endRegion = gene.getStart() - 1;

                    startRegion = Math.max(0, gene.getStart() - maxReadoutConsidered);

                    if (connectedGene != null)
                    {
                        startRegion = Math.max(connectedGene.getEnd() + 1, startRegion);

                    }

                } else
                {

                    startRegion = gene.getEnd() + 1;

                    endRegion = gene.getEnd() + maxReadoutConsidered;

                    if (connectedGene != null)
                    {
                        endRegion = Math.min(connectedGene.getStart() - 1, endRegion);

                    }

                }

                if (endRegion - startRegion + 1 >= windowSize && (gene.getDownstreamGene() != null || gene.isLastInChromosome()))
                {
                    wr.write(gene.getId() + "\t" + gene.getChr() + "\t" + gene.getStrand() + "\t" + geneExpression);

                    if (gene.getExonReadCount() > 0)
                    {

                        readoutEnd = this.getReadoutEnd(reader, gene.getChr(), startRegion, endRegion, gene.getStrand(), strandReversed, geneExpression, threshold, windowSize, useGeneExpression);
                        wr.write("\t" + startRegion + "\t" + endRegion + "\t" + readoutEnd + "\n");

                        if (gene.getStrand() < 0)
                        {
                            wrBed.write(gene.getChr() + "\t" + (readoutEnd - 1) + "\t" + endRegion + "\t" + gene.getSymbol() + "\t1\t-\n");
                        } else
                        {
                            wrBed.write(gene.getChr() + "\t" + (startRegion - 1) + "\t" + readoutEnd + "\t" + gene.getSymbol() + "\t1\t+\n");
                        }

                    } else
                    {
                        wr.write("\tNA\tNA\tNA\n");
                    }
                }

            }
            wr.close();
            wrBed.close();
        } catch (IOException ex)
        {
            Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
        } finally
        {
            try
            {
                wr.close();
            } catch (IOException ex)
            {
                Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    public void getReadoutLevelsInWindows(String bamFile, String outfile, int maxReadoutConsidered, boolean strandReversed, int windowSize, int minOverlap, int minDistanceToDownsteamGene)
    {
        BufferedWriter wr = null;
        BufferedWriter wrBed = null;
        System.out.println(bamFile);
        System.out.println(outfile);
        try
        {
            SamReaderFactory srf = SamReaderFactory.make();
            srf.validationStringency(ValidationStringency.SILENT);
            srf.enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
            SamReader reader = srf.open(new File(bamFile));
            // get reads overlapping this region
            Gene gene;
            int startRegion;
            int endRegion;
            Gene connectedGene;
            int readoutEnd;
            wr = new BufferedWriter(new FileWriter(outfile));
           
            wr.write("Ensembl_ID\tChr\tStrand\tGene_expression\tRegion_start\tRegion_end\twindowsize");

            int maxCount=0;
            
            for(int i=windowSize; i<=maxReadoutConsidered;i+=windowSize)
            {
                wr.write("\t"+i);
                maxCount++;
            }
            wr.write("\n");
            

            double geneExpression;
            double readout;
            int[] counts;
            boolean verbose;

            for (String g : genes.keySet())
            {
                gene = this.genes.get(g);

                // quantify readin;
                geneExpression = (double) gene.getExonReadCount() / gene.getExonLength();

                // quantify readout;
                connectedGene = gene.getDownstreamGene();

                readoutEnd = -1;
                if (gene.getStrand() < 0)
                {
                    endRegion = gene.getStart() - 1;

                    startRegion = Math.max(0, gene.getStart() - maxReadoutConsidered);

                } else
                {

                    startRegion = gene.getEnd() + 1;

                    endRegion = gene.getEnd() + maxReadoutConsidered;

                }

               

                if (gene.getExonReadCount() > 0 && (distances==null || (distances.containsKey(gene.getId()) && distances.get(gene.getId())>=minDistanceToDownsteamGene)))
                {
                    wr.write(gene.getId() + "\t" + gene.getChr() + "\t" + gene.getStrand() + "\t" + geneExpression + "\t"+startRegion+"\t"+endRegion+"\t"+windowSize);
                   
                    verbose=false;
                    counts = this.getReadcountsInWindows(reader, gene.getChr(), startRegion, endRegion, gene.getStrand(), strandReversed, windowSize, minOverlap, verbose);

                    for (int i = 0; i < counts.length; i++)
                    {
                        wr.write("\t" + counts[i]);
                    }
                    if(counts.length<maxCount)
                    {
                        for(int i=counts.length;i<maxCount;i++)
                        {
                             wr.write("\t0");
                        }
                    }
                    
                    wr.write("\n");
                }

                

            }
            wr.close();
          
        } catch (IOException ex)
        {
            Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
        } finally
        {
            try
            {
                wr.close();
            } catch (IOException ex)
            {
                Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    public int getReadoutEnd(SamReader reader, String chr, int start, int end, int strand, boolean strandReversed, double geneExpression, double threshold, int windowSize, boolean useGeneExpression)
    {
        SAMRecord samRecord;
        List<AlignmentBlock> blocks;
        
       
        SAMRecordIterator iterator = reader.query(chr, start, end, false);

        int lengthContained;
        int readstrand;
        HashSet<String> counted = new HashSet<>();

        int[] counts = new int[end - start + 1];

        int startPos;
        int endPos;
        while (iterator.hasNext())
        {
            samRecord = iterator.next();

            readstrand = 1;
            if (samRecord.getReadNegativeStrandFlag())
            {
                if (!samRecord.getReadPairedFlag() || samRecord.getFirstOfPairFlag())
                {
                    if (!strandReversed)
                    {
                        readstrand = -1;
                    }
                } else
                {
                    if (strandReversed)
                    {
                        readstrand = -1;
                    }
                }
            } else
            {
                if (!samRecord.getReadPairedFlag() || samRecord.getFirstOfPairFlag())
                {
                    if (strandReversed)
                    {
                        readstrand = -1;
                    }
                } else
                {
                    if (!strandReversed)
                    {
                        readstrand = -1;
                    }
                }
            }

            if (strand == readstrand)
            {
                blocks = samRecord.getAlignmentBlocks();
                lengthContained = 0;

                for (AlignmentBlock block : blocks)
                {

                    startPos = Math.max(start, block.getReferenceStart());
                    endPos = Math.min(end, block.getReferenceStart() + block.getLength() - 1);

                    for (int i = startPos; i <= endPos; i++)
                    {
                        counts[i - start]++;
                    }

                }

            }

        }
        iterator.close();
        int sum;

        double baseExpression;
        if (strand > 0)
        {
            sum = 0;

            startPos = start;
            endPos = start + windowSize - 1;

            for (int i = startPos; i <= endPos; i++)
            {
                sum += counts[i - start];
            }

            baseExpression = geneExpression;
            if (!useGeneExpression)
            {
                baseExpression = (double) sum / windowSize;
            }

            while ((double) sum / windowSize >= threshold * baseExpression && endPos < end && sum > 0)
            {
                sum -= counts[startPos - start];
                startPos++;
                endPos++;
                sum += counts[endPos - start];
            }

            return endPos;

        } else
        {

            sum = 0;

            startPos = end - windowSize + 1;
            endPos = end;

            for (int i = startPos; i <= endPos; i++)
            {
                sum += counts[i - start];
            }

            baseExpression = geneExpression;
            if (!useGeneExpression)
            {
                baseExpression = (double) sum / windowSize;
            }

            while ((double) sum / windowSize >= threshold * baseExpression && startPos > start && sum > 0)
            {
                sum -= counts[endPos - start];
                startPos--;
                endPos--;
                sum += counts[startPos - start];
            }

            return startPos;

        }

    }

    public int[] getReadcountsInWindows(SamReader reader, String chr, int start, int end, int strand, boolean strandReversed, int windowSize, int minOverlap, boolean verbose)
    {
        SAMRecord samRecord;
        List<AlignmentBlock> blocks;
        

        int nrWindows = (end - start + 1) / windowSize;

        if ((end - start + 1) % windowSize > 0)
        {
            nrWindows++;
        }

        int[] counts = new int[nrWindows];
        int[] starts = new int[nrWindows];
        int[] ends = new int[nrWindows];
        int k;
        
        if ((end - start + 1) % windowSize == 0 || strand > 0)
        {
            for (int i = 0; i < nrWindows; i++)
            {

                starts[i] = start + i * windowSize;
                ends[i] = start + (i + 1) * windowSize - 1;

            }
        }
        else{
            
            for (int i = nrWindows-1; i >=0; i--)
            {
                k=nrWindows-1-i;
                
                starts[i] = Math.max(end - (k+1) * windowSize+1,0);
                ends[i] = Math.max(end- k*windowSize,0);
                
                 if(verbose)
                 {
                     System.out.println(starts[i]+"\t"+ends[i]);
                 }

            }
        }
        
       

        int lengthContained;
        int readstrand;
        HashSet<String> counted;
        SAMRecordIterator iterator;

        for (int i = 0; i < nrWindows; i++)
        {
            counted = new HashSet<>();
            iterator = reader.query(chr, starts[i], ends[i], false);
            

            while (iterator.hasNext())
            {
                samRecord = iterator.next();

                readstrand = 1;
                if (samRecord.getReadNegativeStrandFlag())
                {
                    if (!samRecord.getReadPairedFlag() ||samRecord.getFirstOfPairFlag())
                    {
                        if (!strandReversed)
                        {
                            readstrand = -1;
                        }
                    } else
                    {
                        if (strandReversed)
                        {
                            readstrand = -1;
                        }
                    }
                } else
                {
                    if (!samRecord.getReadPairedFlag() ||samRecord.getFirstOfPairFlag())
                    {
                        if (strandReversed)
                        {
                            readstrand = -1;
                        }
                    } else
                    {
                        if (!strandReversed)
                        {
                            readstrand = -1;
                        }
                    }
                }

                if (strand == readstrand)
                {
                    blocks = samRecord.getAlignmentBlocks();
                    lengthContained = 0;

                    for (AlignmentBlock block : blocks)
                    {
                        lengthContained += Math.max(0, Math.min(ends[i], block.getReferenceStart() + block.getLength() - 1) - Math.max(starts[i], block.getReferenceStart()) + 1);

                    }

                    if (lengthContained >= minOverlap && !counted.contains(samRecord.getReadName()))
                    {
                        counts[i]++;
                        counted.add(samRecord.getReadName());
                        
                        if (verbose && i==nrWindows - 1)
                        {
                            System.out.println(samRecord.toString()+"\t"+lengthContained);
                        }
                    }

                }
            }
            iterator.close();
            
            
        }
        
        if(strand>0)
        {
            return counts;
        }
        else{
            int[] tmp= new int[counts.length];
            
            for(int i=0; i<counts.length;i++)
            {
                tmp[i]=counts[counts.length-1-i];
            }
            return(tmp);
        }

    }

    public int countReads(SamReader reader, String chr, int start, int end, int minOverlap, int strand, boolean strandReversed, boolean strandSpecific, boolean verbose)
    {
        SAMRecord samRecord;
        List<AlignmentBlock> blocks;
        int count = 0;
        
        SAMRecordIterator iterator = reader.query(chr, start, end, false);

        int lengthContained;
        int readstrand;
        HashSet<String> counted = new HashSet<>();
        while (iterator.hasNext())
        {
            samRecord = iterator.next();

            readstrand = 1;
            if (samRecord.getReadNegativeStrandFlag())
            {
                if (!samRecord.getReadPairedFlag() ||samRecord.getFirstOfPairFlag())
                {
                    if (!strandReversed)
                    {
                        readstrand = -1;
                    }
                } else
                {
                    if (strandReversed)
                    {
                        readstrand = -1;
                    }
                }
            } else
            {
                if (!samRecord.getReadPairedFlag() ||samRecord.getFirstOfPairFlag())
                {
                    if (strandReversed)
                    {
                        readstrand = -1;
                    }
                } else
                {
                    if (!strandReversed)
                    {
                        readstrand = -1;
                    }
                }
            }

            if (strand == readstrand || !strandSpecific)
            {
                blocks = samRecord.getAlignmentBlocks();
                lengthContained = 0;

                for (AlignmentBlock block : blocks)
                {
                    lengthContained += Math.max(0, Math.min(end, block.getReferenceStart() + block.getLength() - 1) - Math.max(start, block.getReferenceStart()) + 1);

                   
                }

                if (lengthContained >= minOverlap && !counted.contains(samRecord.getReadName()))
                {
                    count++;
                    counted.add(samRecord.getReadName());
                    if(verbose)
                    {
                        System.out.println(samRecord.toString());
                    }
                }
            }

        }
        iterator.close();
        return count;
    }

    public void readConnections(String... connectionsFiles)
    {
        HashMap<String, Gene> genesWithConnections = new HashMap<String, Gene>();
        try
        {
            for (String connectionsFile : connectionsFiles)
            {
                BufferedReader rd = new BufferedReader(new FileReader(connectionsFile));
                String inline = rd.readLine();
                String[] content;

                String[] info;
                String id;
                String before_id;
                String after_id;
                Gene gene;
                Gene connectedGene;
                while (inline != null)
                {
                    content = inline.split("\t");

                    info = content[8].split(";");

                    before_id = null;
                    after_id = null;
                    id = null;

                    for (String s : info)
                    {
                        if (s.startsWith("ID="))
                        {
                            id = s.substring(3);
                        }
                        if (s.startsWith("before_ID="))
                        {
                            before_id = s.substring(10);
                        }
                        if (s.startsWith("after_ID="))
                        {
                            after_id = s.substring(9);
                        }
                    }

                    gene = genes.get(id);

                    if (!genesWithConnections.containsKey(gene.getId()))
                    {
                        genesWithConnections.put(gene.getId(), gene);
                    }

                    if (content[2].equals("UP"))
                    {
                        if (before_id != null)
                        {

                            if (!before_id.startsWith("NA"))
                            {
                                connectedGene = genes.get(before_id);
                                gene.setUpstreamGene(connectedGene);

                                if (!genesWithConnections.containsKey(connectedGene.getId()))
                                {
                                    genesWithConnections.put(connectedGene.getId(), connectedGene);
                                }
                            } else
                            {
                                gene.setFirstInChromosome(true);
                            }
                        }

                    }
                    if (content[2].equals("DOWN"))
                    {
                        if (after_id != null)
                        {
                            if (!after_id.startsWith("NA"))
                            {
                                connectedGene = genes.get(after_id);
                                gene.setDownstreamGene(connectedGene);

                                if (!genesWithConnections.containsKey(connectedGene.getId()))
                                {
                                    genesWithConnections.put(connectedGene.getId(), connectedGene);
                                }
                            } else
                            {
                                gene.setLastInChromosome(true);
                            }
                        }

                    }

                    inline = rd.readLine();
                }
                rd.close();
            }

        } catch (IOException ex)
        {
            Logger.getLogger(ReadOutFromBamCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }

        System.out.println(genesWithConnections.size());
        this.genes = genesWithConnections;
    }

    
}
