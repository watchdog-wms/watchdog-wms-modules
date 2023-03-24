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
import java.util.HashSet;

/**
 *
 * @author friedel
 */
public class GeneInformation implements Comparable<GeneInformation>
{
    private String geneId;
    private String geneName;
    private String chr;
    private int strand;
    private int geneStart;
    private int geneEnd;
    private int last3PrimeUTRStart=-1;
    private int last3PrimeUTREnd=-1;
    private boolean hasUTR;
    private ArrayList<Exon> UTRs;
    
    private GeneInformation geneLink;
  

    public GeneInformation()
    {
        hasUTR=false;
        //this.geneLink=null;
    }

    public GeneInformation(String geneId,String geneName,  String chr, int strand, int geneStart, int geneEnd)
    {
        this();
        this.geneId = geneId;
        this.geneName=geneName;
        this.chr = chr;
        this.strand = strand;
        this.geneStart = geneStart;
        this.geneEnd = geneEnd;
        this.UTRs= new ArrayList<>();
     
        
    }
    
   
    
    public void add3PrimeUTR(int start, int end)
    {
        hasUTR=true;
        Exon exon= new Exon();
        exon.setStart(start);
        exon.setEnd(end);
        this.UTRs.add(exon);
        if(last3PrimeUTRStart <0)
        {
            this.last3PrimeUTRStart=start;
            this.last3PrimeUTREnd=end;
        }
        else{
            
            if(strand >0)
            {
                if(end > this.last3PrimeUTREnd || (end==this.last3PrimeUTREnd && start < this.last3PrimeUTRStart))
                {
                    this.last3PrimeUTRStart = start;
                    this.last3PrimeUTREnd = end;
                   
                }
                
            }
            else{
                if(start < this.last3PrimeUTRStart || (start == this.last3PrimeUTRStart && end > this.last3PrimeUTREnd))
                {
                    this.last3PrimeUTRStart = start;
                    this.last3PrimeUTREnd = end;
                     
                }
            }
        }
    }
    
    public String getUTRSegments()
    {
        Exon[] exons= this.UTRs.toArray(new Exon[0]);
        Arrays.sort(exons);
       
        HashSet<Integer> cuts= new HashSet<>();
        StringBuilder sb= new StringBuilder();
        int start;
        int end;
        String strandinfo="+";
        if(this.strand<0)
        {
             strandinfo="-";
        }
        for(int i=0; i< exons.length;i++)
        {
            cuts.add(exons[i].getStart());
            cuts.add(exons[i].getEnd()+1);
            if(i== exons.length-1 || exons[i].getEnd() <exons[i+1].getStart())
            {
                
                Integer[] pos= cuts.toArray(new Integer[0]);
                Arrays.sort(pos);
                for(int j=0; j< pos.length-1;j++)
                {
                    start=pos[j];
                    end=pos[j+1];
                    
                    end--;
                    
                    sb.append("chr" + this.getChr() + "\tPRED\tUTRSEG\t" + start + "\t" + end + "\t.\t"+strandinfo+"\t.\tID=" + this.geneId +"_"+start+"_"+end+"\n");

                }
              
                cuts= new HashSet<>();
            }
        }
        return sb.toString();
    }
    
    public String getKey()
    {
        return this.geneId;
    }

    @Override
    public String toString()
    {
        return this.geneId+"\t"+this.geneName+"\t"+ this.chr+"\t"+ this.strand+"\t"+ this.geneStart+"\t"+this.geneEnd+"\t"+this.last3PrimeUTRStart+"\t"+this.last3PrimeUTREnd;
    }

    @Override
    public int compareTo(GeneInformation o)
    {
        if(!this.chr.equals(o.chr))
        {
           // return this.chr.compareTo(o.chr);
            return this.compareChromosomes(this.chr, o.chr);
        }
        else
        {
            if(this.geneStart < o.geneStart)
            {
                return -1;
            }
            else if(this.geneStart > o.geneStart)
            {
                return 1;
            }
            else{
                return this.geneEnd-o.geneEnd;
            }
        }
    }
    
    public int compareChromosomes(String chr1, String chr2)
    {
        if(chr1.equals(chr2))
        {
            return 0;
        }
        else{
            try
            {
                
                int chr1ToInt = Integer.parseInt(chr1);
                int chr2ToInt = Integer.parseInt(chr2);
                
                return chr1ToInt-chr2ToInt;
            }
            catch(NumberFormatException exp)
            {
                return chr1.compareTo(chr2);
            }
        }
    }

    public String getChr()
    {
        return chr;
    }

    public int getStrand()
    {
        return strand;
    }

    public int getLast3PrimeUTREnd()
    {
        return last3PrimeUTREnd;
    }

    public int getLast3PrimeUTRStart()
    {
        return last3PrimeUTRStart;
    }

    public int getGeneEnd()
    {
        return geneEnd;
    }

    public int getGeneStart()
    {
        return geneStart;
    }

    public String getGeneId()
    {
        return geneId;
    }

    public void setGeneName(String geneName)
    {
        this.geneName = geneName;
    }

    public String getGeneName()
    {
        return geneName;
    }

    public boolean hasUTR()
    {
        return hasUTR;
    }

    public void setChr(String chr)
    {
        this.chr = chr;
    }

    public void setGeneEnd(int geneEnd)
    {
        this.geneEnd = geneEnd;
    }

    public void setGeneId(String geneId)
    {
        this.geneId = geneId;
    }

    public void setGeneStart(int geneStart)
    {
        this.geneStart = geneStart;
    }

    public void setStrand(int strand)
    {
        this.strand = strand;
    }
    
    
    public void setGeneLink(GeneInformation gene)
    {
        this.geneLink=gene;
    }
    
    public GeneInformation getGeneLink()
    {
        return this.geneLink;
    }
    
    
    public boolean hasGeneLink()
    {
        return this.geneLink != null;
    }
    
}
