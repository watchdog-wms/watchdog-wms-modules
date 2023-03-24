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


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;

/**
 *
 * @author friedel
 */
public class PeakSummary
{

    private HashMap<String, ArrayList<Peak>> peaks;
    private Peak promoter;

    private static final String[] PROMOTERKEYS = new String[]
    {
        "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)"
    };

    private Gene gene;

    private ArrayList<Peak> downstreamPeaks = new ArrayList<>();
    
    private ArrayList<Peak> genePeaks = new ArrayList<>();
    
    private int downstreamEnd;
    private int totalLengthCovered;
    private double weightedLengthCovered;
    private double minScore;
    private double maxScore;
    private int maxLength;
    private int minLength;
    
    private int inGeneLength;

    public PeakSummary()
    {
        this.peaks = new HashMap<>();
    }

    public void setGene(Gene gene)
    {
        this.gene = gene;
    }

    public Gene getGene()
    {
        return gene;
    }

    public Peak getPromoter()
    {
        return promoter;
    }

    public void addPeak(Peak peak)
    {
        if (!this.peaks.containsKey(peak.getAnnotation()))
        {
            this.peaks.put(peak.getAnnotation(), new ArrayList<Peak>());

        }
        this.peaks.get(peak.getAnnotation()).add(peak);
    }

    public void summarize(int maxDistance)
    {
        summarize( maxDistance, false);
    }
    public void summarize(int maxDistance, boolean countOverlapping)
    {
        // find promoter;
        promoter = null;
        int i = 0;
        ArrayList<Peak> cands;
        while (promoter == null && i < PROMOTERKEYS.length)
        {

            cands = peaks.get(PROMOTERKEYS[i]);
            if (cands != null)
            {
                for (Peak c : cands)
                {
                    if (promoter == null)
                    {
                        promoter = c;
                    } else
                    {
                        if (Math.abs(promoter.getDistanceToTSS()) > Math.abs(c.getDistanceToTSS()))
                        {
                            promoter = c;
                        } else if (promoter.getDistanceToTSS() == c.getDistanceToTSS())
                        {
                            if (promoter.getScore() < c.getScore())
                            {
                                promoter = c;
                            }
                        }
                    }
                }
            }

            i++;
        }

        // get downstream region
        int strand;

        downstreamEnd = 0;
        totalLengthCovered = 0;
        weightedLengthCovered = 0;
        minScore = Double.MAX_VALUE;
        maxScore = 0;
        
        this.maxLength = 0;
        this.minLength=Integer.MAX_VALUE;
        
        this.downstreamPeaks = new ArrayList<Peak>();
        boolean add;
        int lengthDownStream;
        if (gene != null)
        {
            strand = gene.getStrand();
            for (String anno : this.peaks.keySet())
            {

                cands = peaks.get(anno);
                if (cands != null)
                {
                    for (Peak c : cands)
                    {
                        add = false;
                        lengthDownStream = 0;
                        if(!countOverlapping)
                        {
                            if (strand > 0 && c.getEnd() > gene.getEnd() && c.getStart() > gene.getEnd())
                            {
                                add = true;

                            } else if (strand < 0 && c.getStart() < gene.getStart() && c.getEnd() < gene.getStart())
                            {
                                add = true;

                            }
                        }
                        else{
                            if (strand > 0 && c.getEnd() > gene.getEnd())
                            {
                                add = true;

                            } else if (strand < 0 && c.getStart() < gene.getStart())
                            {
                                add = true;

                            }
                        }
                            
                        if (add)
                        {
                            this.downstreamPeaks.add(c);
                            

                        }
                    }
                }
            }
            
            Peak[] downpeaks=this.downstreamPeaks.toArray(new Peak[0]);
           
         
            
            this.downstreamPeaks = new ArrayList<Peak>();

            int startDownstream=0;
            int endDownstream=0;
            HashMap<Integer, Double> downStreamPos= new HashMap<Integer, Double>();
            
            int lastEnd=gene.getEnd();
            downstreamEnd=gene.getEnd();
            if(strand<0)
            {
                lastEnd=gene.getStart();
                downstreamEnd=gene.getStart();
            }
            
            
            for (Peak c : downpeaks)
            {
                
                add = false;
                lengthDownStream = 0;
                if (strand > 0 && c.getEnd() > gene.getEnd())
                {
                    lastEnd=c.getEnd();
                    add = true;
                    lengthDownStream = c.getEnd() - Math.max(gene.getEnd(), c.getStart() - 1);
                    startDownstream=Math.max(gene.getEnd()+1, c.getStart());
                    endDownstream= c.getEnd(); 
                    if (c.getEnd() > downstreamEnd)
                    {
                        downstreamEnd = c.getEnd();
                    }

                } else if (strand < 0 && c.getStart() < gene.getStart() )
                {
                    lastEnd=c.getStart();
                    add = true;
                    lengthDownStream = Math.min(gene.getStart(), c.getEnd() + 1) - c.getStart();
                    startDownstream=c.getStart();
                    endDownstream=Math.min(gene.getStart()+1, c.getEnd());
                    
                    if (c.getStart() < downstreamEnd)
                    {
                        downstreamEnd = c.getStart();
                    }
                }
                if (add)
                {
                    this.downstreamPeaks.add(c);
                    
                   
                    
                    for(int j=startDownstream;j<=endDownstream;j++)
                    {
                        
                        if(downStreamPos.containsKey(j))
                        {
                            downStreamPos.put(j, Math.max(downStreamPos.get(j), c.getScore()));
                        }
                        else{
                             downStreamPos.put(j, c.getScore());
                        }
                    }

                    if (maxScore < c.getScore())
                    {
                        maxScore = c.getScore();
                    }
                    if (minScore > c.getScore())
                    {
                        minScore = c.getScore();
                    }
                    
                    int clength= c.getEnd()-c.getStart()+1;
                    if (maxLength < clength)
                    {
                        maxLength = clength;
                    }
                    if (minLength > clength)
                    {
                        minLength = clength;
                    }
                    

                }
            }
            
            for(Integer pos:downStreamPos.keySet())
            {
                
                totalLengthCovered++;
                weightedLengthCovered+=downStreamPos.get(pos);
            }

        }
    }
    
    public void summarizeInGenes()
    {
        ArrayList<Peak> cands;
        int strand;
        this.genePeaks= new ArrayList<Peak>();
        if (gene != null)
        {
            strand = gene.getStrand();
            for (String anno : this.peaks.keySet())
            {

                cands = peaks.get(anno);
                if (cands != null)
                {
                    for (Peak c : cands)
                    {
                         this.genePeaks.add(c);
                         this.inGeneLength+=Math.min(c.getEnd(), this.gene.getEnd())-Math.max(c.getStart(), this.gene.getStart())+1;
                    }
                }
            }
        }
    }
    
   

    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        sb.append(this.gene.getId()).append("\t");
        sb.append(this.gene.getSymbol()).append("\t");
        sb.append(this.gene.getChr()).append("\t");
        sb.append(this.gene.getStrand()).append("\t");
        sb.append(this.gene.getStart()).append("\t");
        sb.append(this.gene.getEnd()).append("\t");
        sb.append(this.promoterToString()).append("\t");
        sb.append(this.downstreamEnd).append("\t");
        sb.append(this.maxScore).append("\t");
        sb.append(this.minScore).append("\t");
        sb.append(this.maxLength).append("\t");
        sb.append(this.minLength).append("\t");
        sb.append(this.totalLengthCovered).append("\t");
        sb.append(this.weightedLengthCovered).append("\t");
        sb.append(this.downstreamToString());
        return sb.toString();
    }
    
    public String toStringInGenes()
    {
        StringBuilder sb = new StringBuilder();
        sb.append(this.gene.getId()).append("\t");
        sb.append(this.gene.getSymbol()).append("\t");
        sb.append(this.gene.getChr()).append("\t");
        sb.append(this.gene.getStrand()).append("\t");
        sb.append(this.gene.getStart()).append("\t");
        sb.append(this.gene.getEnd()).append("\t");
        sb.append(this.inGeneLength).append("\t");
        sb.append((double)this.inGeneLength/(this.gene.getEnd()-this.gene.getStart()+1)).append("\t");
        sb.append(this.inGeneToString());
        return sb.toString();
    }
    
    
    
    
    public String getBEDDownstream()
    {
        StringBuilder sb = new StringBuilder();

        Peak[] downpeaks=this.downstreamPeaks.toArray(new Peak[0]);
         
        Arrays.sort(downpeaks, new ComparatorPositiveStrand());
            
       
        for (Peak c : downpeaks)
        {
            sb.append(c.getChr()).append("\t").append(c.getStart()).append("\t").append(c.getEnd()).append("\t").append(c.getId()).append("_").append(this.getGene().getSymbol()).append("\t").append(c.getScore()).append("\n");
        }
       

        return sb.toString();
    }
    
      public String getBEDInGene()
    {
        StringBuilder sb = new StringBuilder();

        Peak[] genepeaks=this.genePeaks.toArray(new Peak[0]);
         
        Arrays.sort(genepeaks, new ComparatorPositiveStrand());
            
       
        for (Peak c : genepeaks)
        {
            sb.append(c.getChr()).append("\t").append(c.getStart()).append("\t").append(c.getEnd()).append("\t").append(c.getId()).append("_").append(this.getGene().getSymbol()).append("\t").append(c.getScore()).append("\n");
        }
       

        return sb.toString();
    }
    
    
   
    
    public String getDownStreamRegionAsBED()
    {
        StringBuilder sb = new StringBuilder();
        int strand;
        String chr;
        int downstreamLength;
    
        if (gene != null)
        {
            strand = gene.getStrand();
            chr=gene.getChr();
            sb.append(chr).append("\t");
            if (strand > 0)
            {
                downstreamLength=this.downstreamEnd-gene.getEnd();
                sb.append(gene.getEnd()).append("\t");
                sb.append(this.downstreamEnd).append("\t");
                
            }
            else{
                
                downstreamLength=gene.getStart()-this.downstreamEnd;
                sb.append(this.downstreamEnd).append("\t");
                sb.append(gene.getStart()).append("\t");
               
            }
            sb.append(gene.getSymbol()).append("\t");
            sb.append(this.totalLengthCovered).append("\t");
            
            if(strand > 0)
            {
                sb.append("+\n");
            }
            else{
                sb.append("-\n");
            }
            
        }
        return sb.toString();
    }         
        
    public String promoterToString()
    {
        if (promoter != null)
        {
            return this.promoter.getId() + "\t" + this.promoter.getChr() + "\t" + this.promoter.getStart()
                    + "\t" + this.promoter.getEnd() + "\t" + this.promoter.getScore() +"\t"+ (this.promoter.getEnd()-this.promoter.getStart()+1)+"\t" + this.promoter.getAnnotation() + "\t" + this.promoter.getDistanceToTSS();
        } else
        {
            return "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
        }
    }
    
    

    public String downstreamToString()
    {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < this.downstreamPeaks.size(); i++)
        {
            sb.append(this.downstreamPeaks.get(i).getId());
            if (i < this.downstreamPeaks.size() - 1)
            {
                sb.append(",");
            }
        }

        if (this.downstreamPeaks.isEmpty())
        {
            sb.append("NA");
        }

        return sb.toString();
    }
    
     public String inGeneToString()
    {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < this.genePeaks.size(); i++)
        {
            sb.append(this.genePeaks.get(i).getId());
            if (i < this.genePeaks.size() - 1)
            {
                sb.append(",");
            }
        }

        if (this.genePeaks.isEmpty())
        {
            sb.append("NA");
        }

        return sb.toString();
    }
    
    public static class ComparatorPositiveStrand implements Comparator<Peak>
    {

        @Override
        public int compare(Peak o1, Peak o2)
        {
            if(o1.getChr().equals(o2.getChr()))
            {
                return o1.getStart()-o2.getStart();
            }
            else{
                return o1.getChr().compareTo(o2.getChr());
            }
        }
        
    }
    
    public static class ComparatorNegativeStrand implements Comparator<Peak>
    {

        @Override
        public int compare(Peak o1, Peak o2)
        {
            if(o1.getChr().equals(o2.getChr()))
            {
                return o2.getEnd()-o1.getEnd();
            }
            else{
                return o1.getChr().compareTo(o2.getChr());
            }
        }
        
    }
    
    public boolean isEmpty()
    {
        return this.peaks.isEmpty();
    }

    public HashMap<String, ArrayList<Peak>> getPeaks()
    {
        return peaks;
    }

    

    
}
