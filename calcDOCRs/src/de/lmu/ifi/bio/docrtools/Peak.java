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
public class Peak implements Comparable<Peak>
{
    private String chr;
    private int start;
    private int end;
    
    private double score;
    private String id;
    private String annotation;
    private int distanceToTSS;

    public Peak(String chr, int start, int end, double score, String id, String annotation, int distanceToTSS)
    {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.score = score;
        this.id = id;
        this.annotation = annotation;
        this.distanceToTSS = distanceToTSS;
    }

    public String getAnnotation()
    {
        return annotation;
    }

    public String getChr()
    {
        return chr;
    }

    public int getDistanceToTSS()
    {
        return distanceToTSS;
    }

    public int getEnd()
    {
        return end;
    }


    public String getId()
    {
        return id;
    }

    public double getScore()
    {
        return score;
    }

    public int getStart()
    {
        return start;
    }

    @Override
    public String toString()
    {
        return this.getId()+"\t"+this.getAnnotation()+"\t"+this.getChr()+"\t"+this.getStart()+"\t"+this.getEnd();
    }

    @Override
    public int compareTo(Peak o)
    {
        
        if (this.getChr().equals(o.getChr()))
        {
            return this.getStart() - o.getStart();
        } else
        {
            return this.getChr().compareTo(o.getChr());
        }
        
    }

    

    
    
    
}
