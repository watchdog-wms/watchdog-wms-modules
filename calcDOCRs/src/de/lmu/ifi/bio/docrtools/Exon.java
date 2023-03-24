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

/**
 *
 * @author friedel
 */
public class Exon implements Comparable<Exon>
{
    private String geneid;
    private String exonid;
    private int strand;
    private String chr;
    private int start;
    private int end;
    private String sequence;
    

    public Exon()
    {
    }

    
    public Exon(String geneid, String exonid, int strand, String chr, int start, int end, String sequence)
    {
        this.geneid = geneid;
        this.exonid = exonid;
        this.strand = strand;
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.sequence = sequence;
    }

    public String getExonid()
    {
        return exonid;
    }

    public String getGeneid()
    {
        return geneid;
    }

    public String getChr()
    {
        return chr;
    }

    public int getEnd()
    {
        return end;
    }

    public String getSequence()
    {
        return sequence;
    }

    public int getStart()
    {
        return start;
    }

    public int getStrand()
    {
        return strand;
    }
    
    
    public int getLength()
    {
        return this.sequence.length();
    }

    public void setStart(int start)
    {
        this.start = start;
    }

    public void setEnd(int end)
    {
        this.end = end;
    }

    @Override
    public int compareTo(Exon o)
    {
        if(this.start!= o.start)
        {
            return this.start-o.start;
        }
        else{
            return this.end-o.end;
        }
    }

    @Override
    public String toString()
    {
        return this.chr+"\t"+this.strand+"\t"+this.start+"\t"+this.end;
    }
    
    
    
}
