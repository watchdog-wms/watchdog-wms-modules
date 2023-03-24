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
public class Gene
{
    private String id;
    private String symbol;
    private int strand;
    private String chr;
    private int start;
    private int end;
    
    private Gene upstreamGene;
    private Gene downstreamGene;
    
    private int exonLength;
    private int exonReadCount;
    
    private boolean firstInChromosome;
    private boolean lastInChromosome;

    public Gene(String id, String symbol, int strand, String chr, int start, int end)
    {
        this.id = id;
        this.symbol = symbol;
        this.strand = strand;
        this.chr = chr;
        this.start = start;
        this.end = end;

    }

    public String getChr()
    {
        return chr;
    }

    public int getEnd()
    {
        return end;
    }

    public String getId()
    {
        return id;
    }

    public int getStart()
    {
        return start;
    }

    public Gene getDownstreamGene()
    {
        return downstreamGene;
    }

    public Gene getUpstreamGene()
    {
        return upstreamGene;
    }

    public int getStrand()
    {
        return strand;
    }

    public String getSymbol()
    {
        return symbol;
    }

    public void setChr(String chr)
    {
        this.chr = chr;
    }

    public void setDownstreamGene(Gene downstreamGene)
    {
        this.downstreamGene = downstreamGene;
    }

    public void setEnd(int end)
    {
        this.end = end;
    }

    public void setId(String id)
    {
        this.id = id;
    }

    public void setStart(int start)
    {
        this.start = start;
    }

    public void setStrand(int strand)
    {
        this.strand = strand;
    }

    public void setSymbol(String symbol)
    {
        this.symbol = symbol;
    }

    public void setUpstreamGene(Gene upstreamGene)
    {
        this.upstreamGene = upstreamGene;
    }

    public void setExonLength(int exonLength)
    {
        this.exonLength = exonLength;
    }

    public void setExonReadCount(int exonReadCount)
    {
        this.exonReadCount = exonReadCount;
    }

    public int getExonLength()
    {
        return exonLength;
    }

    public int getExonReadCount()
    {
        return exonReadCount;
    }
    
    
    public void reInitialize()
    {
       this.setExonReadCount(0);
    }

    public void setFirstInChromosome(boolean firstInChromosome)
    {
        this.firstInChromosome = firstInChromosome;
    }

    public void setLastInChromosome(boolean lastInChromosome)
    {
        this.lastInChromosome = lastInChromosome;
    }

    public boolean isFirstInChromosome()
    {
        return firstInChromosome;
    }

    public boolean isLastInChromosome()
    {
        return lastInChromosome;
    }
    
    
    
    
}
