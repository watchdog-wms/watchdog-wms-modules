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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author friedel
 */
public class UTRCalculator
{

    private HashMap<String, GeneInformation> genes;

    public UTRCalculator()
    {
        this.genes = new HashMap<>();
    }

    public HashMap<String, GeneInformation> getGenes()
    {
        return genes;
    }

    public void readGeneInformation(String infile)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String inline = rd.readLine();
            String[] content;
            GeneInformation gene;
            while (inline != null)
            {
                if (!inline.startsWith("Ensembl"))
                {
                    content = inline.split("\t");

                    gene = new GeneInformation(content[0], content[8], content[2], Integer.parseInt(content[7]), Integer.parseInt(content[3]), Integer.parseInt(content[4]));
                    if (this.genes.containsKey(gene.getKey()))
                    {
                        gene = this.genes.get(gene.getKey());
                    } else
                    {
                        this.genes.put(gene.getKey(), gene);
                    }

                    if (content.length > 11)
                    {
                        if (!content[11].equals("") && !content[12].equals(""))
                        {
                            gene.add3PrimeUTR(Integer.parseInt(content[11]), Integer.parseInt(content[12]));
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

    public void readGeneInformation(String infile, int strand, HashSet<String> exclude)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String inline = rd.readLine();
            String[] content;
            GeneInformation gene;
            String biotype;

            while (inline != null)
            {
                if (!inline.startsWith("Ensembl"))
                {
                    content = inline.split("\t");

                    gene = new GeneInformation(content[0], content[8], content[2], Integer.parseInt(content[7]), Integer.parseInt(content[3]), Integer.parseInt(content[4]));
                    biotype = content[content.length - 1];

                    if (gene.getStrand() * strand > 0 && !exclude.contains(biotype))
                    {
                        if (this.genes.containsKey(gene.getKey()))
                        {
                            gene = this.genes.get(gene.getKey());
                        } else
                        {
                            this.genes.put(gene.getKey(), gene);
                        }
                        if (content.length > 11 && !content[11].equals("") && !content[12].equals(""))
                        {
                            gene.add3PrimeUTR(Integer.parseInt(content[11]), Integer.parseInt(content[12]));
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

    public void readGeneInformation(String infile, int strand, HashSet<String> exclude, String patternChromosome)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String inline = rd.readLine();
            String[] content;
            GeneInformation gene;
            String biotype;

            while (inline != null)
            {
                if (!inline.startsWith("Ensembl"))
                {
                    content = inline.split("\t");

                    gene = new GeneInformation(content[0], content[8], content[2], Integer.parseInt(content[7]), Integer.parseInt(content[3]), Integer.parseInt(content[4]));
                    biotype = content[content.length - 1];

                    if (gene.getStrand() * strand > 0 && !exclude.contains(biotype) && gene.getChr().matches(patternChromosome))
                    {
                        if (this.genes.containsKey(gene.getKey()))
                        {
                            gene = this.genes.get(gene.getKey());
                        } else
                        {
                            this.genes.put(gene.getKey(), gene);
                        }
                        if (content.length > 11 && !content[11].equals("") && !content[12].equals(""))
                        {
                            gene.add3PrimeUTR(Integer.parseInt(content[11]), Integer.parseInt(content[12]));
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

    public void readGeneInformationFromGTF(String infile, int strand, HashSet<String> exclude, String patternChromosome)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String inline = rd.readLine();
            String[] content;
            GeneInformation gene;
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
                        chr = chr.replace("chr", "");
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

                        gene = new GeneInformation(geneId, geneName, chr, geneStrand, geneStart, geneEnd);

                        if (gene.getStrand() * strand > 0 && !exclude.contains(biotype) && gene.getChr().matches(patternChromosome))
                        {
                            if (!this.genes.containsKey(gene.getKey()))
                            {
                                this.genes.put(gene.getKey(), gene);
                            }

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

    public void readGeneInformationFromGTF(String infile, HashSet<String> exclude, String patternChromosome)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String inline = rd.readLine();
            String[] content;
            GeneInformation gene;
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
                        chr = chr.replace("chr", "");
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

                        gene = new GeneInformation(geneId, geneName, chr, geneStrand, geneStart, geneEnd);

                       
                        if (!exclude.contains(biotype) && gene.getChr().matches(patternChromosome))
                        {
                            if (!this.genes.containsKey(gene.getKey()))
                            {
                                this.genes.put(gene.getKey(), gene);
                            }

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

    public void assignPTT(GeneInformation[] ptts, int maxDistToGene, String outfile, int strand, int allowedOverlap)
    {
        try
        {
            GeneInformation[] genes_sorted = this.genes.values().toArray(new GeneInformation[0]);
            GeneInformation[] elements = new GeneInformation[ptts.length + genes_sorted.length];
            System.arraycopy(genes_sorted, 0, elements, 0, genes_sorted.length);
            System.arraycopy(ptts, 0, elements, genes_sorted.length, ptts.length);
            Arrays.sort(elements);

            HashSet<GeneInformation> previousGenes = new HashSet<GeneInformation>();
            HashSet<GeneInformation> tmp;

            int minDist;
            GeneInformation bestGene;
            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));

            GeneInformation lastPTT = null;
            GeneInformation ptt;

            GeneInformation link;
            int prevDist;
            for (GeneInformation element : elements)
            {
                if (element.getGeneName().equals("PTT"))
                {
                    ptt = element;

                    if (strand < 0)
                    {
                        if (lastPTT != null)
                        {
                            ptt = lastPTT;
                        } else
                        {
                            ptt = null;
                        }
                    }
                    tmp = new HashSet<GeneInformation>();

                    if (ptt != null)
                    {

                        minDist = Integer.MAX_VALUE;
                        bestGene = null;

                        for (GeneInformation g : previousGenes)
                        {
                            if (ptt.getChr().equals(g.getChr()))
                            {
                                if (strand > 0 && g.getGeneEnd() < ptt.getGeneStart() + allowedOverlap)
                                {
                                    if (ptt.getGeneStart() - g.getGeneEnd() <= Math.min(maxDistToGene, minDist))
                                    {
                                        minDist = ptt.getGeneStart() - g.getGeneEnd();
                                        bestGene = g;
                                    }
                                } else if (strand < 0 && g.getGeneStart() > ptt.getGeneEnd() - allowedOverlap)
                                {
                                    if (g.getGeneStart() - ptt.getGeneEnd() <= Math.min(maxDistToGene, minDist))
                                    {

                                        minDist = g.getGeneStart() - ptt.getGeneEnd();
                                        bestGene = g;
                                    }
                                } else
                                {
                                    tmp.add(g);
                                }
                            }
                        }

                        if (bestGene != null)
                        {
                            wr.write(bestGene.toString() + "\t" + ptt.getChr() + "\t" + ptt.getGeneStart() + "\t" + ptt.getGeneEnd() + "\n");

                            if (ptt.hasGeneLink())
                            {
                                link = ptt.getGeneLink();

                                if (link.getStrand() > 0)
                                {
                                    prevDist = ptt.getGeneStart() - link.getGeneEnd();
                                } else
                                {
                                    prevDist = link.getGeneStart() - ptt.getGeneEnd();
                                }

                                if (prevDist > minDist)
                                {
                                    ptt.setGeneLink(bestGene);
                                }

                            } else
                            {
                                ptt.setGeneLink(bestGene);
                            }

                        } else
                        {
                           
                        }

                    }

                    previousGenes = tmp;
                    lastPTT = element;
                } else
                {
                    previousGenes.add(element);
                }

            }

            if (lastPTT != null && strand < 0)
            {
                ptt = lastPTT;

                minDist = Integer.MAX_VALUE;
                bestGene = null;

                for (GeneInformation g : previousGenes)
                {
                    if (ptt.getChr().equals(g.getChr()))
                    {
                        if (strand < 0 && g.getGeneStart() > ptt.getGeneEnd() - allowedOverlap)
                        {
                            if (g.getGeneStart() - ptt.getGeneEnd() <= Math.min(maxDistToGene, minDist))
                            {

                                minDist = g.getGeneStart() - ptt.getGeneEnd();
                                bestGene = g;
                            }
                        }

                    }
                }

                if (bestGene != null)
                {
                    wr.write(bestGene.toString() + "\t" + ptt.getChr() + "\t" + ptt.getGeneStart() + "\t" + ptt.getGeneEnd() + "\n");

                    if (ptt.hasGeneLink())
                    {
                        link = ptt.getGeneLink();

                        if (link.getStrand() > 0)
                        {
                            prevDist = ptt.getGeneStart() - link.getGeneEnd();
                        } else
                        {
                            prevDist = link.getGeneStart() - ptt.getGeneEnd();
                        }

                        if (prevDist > minDist)
                        {
                            ptt.setGeneLink(bestGene);
                        }

                    } else
                    {
                        ptt.setGeneLink(bestGene);
                    }
                }
            }
            wr.close();
        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void printPTTAssignment(GeneInformation[] ptts, String outfile)
    {
        try
        {
            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));

            int countUnassigned = 0;

            for (GeneInformation ptt : ptts)
            {
                if (ptt.hasGeneLink())
                {
                    wr.write(ptt.getChr() + "\t" + ptt.getGeneStart() + "\t" + ptt.getGeneEnd() + "\t" + ptt.getGeneLink().toString() + "\n");
                } else
                {
                    countUnassigned++;
                    System.out.println(ptt.getChr() + ":" + ptt.getGeneStart() + "-" + ptt.getGeneEnd());
                }
            }

            wr.close();
            System.out.println(countUnassigned);
        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void getConnections(String outfile, String fileDistances, String fileGenePairs, int maxLength, int minLength, int strand)
    {
        try
        {
            GeneInformation[] g = this.genes.values().toArray(new GeneInformation[0]);
            Arrays.sort(g);

            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
            int rightMostEnd = -1;
            String chr = null;
            BufferedWriter wr2 = new BufferedWriter(new FileWriter(fileDistances));
            BufferedWriter wr3 = new BufferedWriter(new FileWriter(fileGenePairs));

            GeneInformation rightMostGene = null;
            GeneInformation gene;

            for (int i = 0; i < g.length; i++)
            {
                

                if (chr != null && (i + 1 < g.length && !g[i].getChr().equals(g[i + 1].getChr())) || (i + 1 == g.length))
                {
                    if (rightMostGene.getGeneEnd() >= g[i].getGeneEnd())
                    {
                        gene = rightMostGene;
                    } else
                    {
                        gene = g[i];
                    }
                    int enddown = gene.getGeneEnd() + maxLength;
                    int startdown = gene.getGeneEnd() + 1;

                    if (strand < 0)
                    {
                        wr.write("chr" + gene.getChr() + "\tPRED\tUP\t" + startdown + "\t" + enddown + "\t.\t-\t.\tID=" + gene.getGeneId() + ";gene_id=" + gene.getGeneName() + ";before_ID=NA" + "\n");

                    } else
                    {
                        wr.write("chr" + gene.getChr() + "\tPRED\tDOWN\t" + startdown + "\t" + enddown + "\t.\t+\t.\tID=" + gene.getGeneId() + ";gene_id=" + gene.getGeneName() + ";after_ID=NA" + "\n");

                    }
                }

                if (chr == null || !chr.equals(g[i].getChr()) || g[i].getGeneStart() > rightMostEnd)
                {

                    if (g[i].getGeneEnd() > rightMostEnd || !chr.equals(g[i].getChr()))
                    {
                        rightMostEnd = g[i].getGeneEnd();
                        rightMostGene = g[i];
                        chr = g[i].getChr();
                    }
                    wr.write(g[i].getUTRSegments());

                    if (strand < 0)
                    {
                        wr.write("chr" + g[i].getChr() + "\tPRED\tANTI\t" + g[i].getGeneStart() + "\t" + g[i].getGeneEnd() + "\t.\t+\t.\tID=" + g[i].getGeneId() + ";gene_id=" + g[i].getGeneName() + "\n");

                    } else
                    {
                        wr.write("chr" + g[i].getChr() + "\tPRED\tANTI\t" + g[i].getGeneStart() + "\t" + g[i].getGeneEnd() + "\t.\t-\t.\tID=" + g[i].getGeneId() + ";gene_id=" + g[i].getGeneName() + "\n");

                    }

                    int next = i + 1;

                    while (next < g.length && g[next].getGeneStart() <= g[i].getGeneEnd() && g[next].getGeneEnd() <= g[i].getGeneEnd() && g[next].getChr().equals(g[i].getChr()))
                    {
                        next++;

                    }

                    if (next < g.length && g[next].getChr().equals(g[i].getChr()))
                    {
                        if (g[next].getGeneStart() > g[i].getGeneEnd() + minLength)
                        {
                            int intervalLength = g[next].getGeneStart() - g[i].getGeneEnd() - 1;
                            int updownsize = Math.min(intervalLength / 3, maxLength);
                            int enddown = g[i].getGeneEnd() + updownsize;
                            int startdown = g[i].getGeneEnd() + 1;

                            int endup = g[next].getGeneStart() - 1;
                            int startup = g[next].getGeneStart() - updownsize;

                            if (strand < 0)
                            {
                                wr.write("chr" + g[i].getChr() + "\tPRED\tUP\t" + startdown + "\t" + enddown + "\t.\t-\t.\tID=" + g[i].getGeneId() + ";gene_id=" + g[i].getGeneName() + ";before_ID=" + g[next].getGeneId() + ";before_gene_id=" + g[next].getGeneName() + "\n");
                                wr.write("chr" + g[next].getChr() + "\tPRED\tDOWN\t" + startup + "\t" + endup + "\t.\t-\t.\tID=" + g[next].getGeneId() + ";gene_id=" + g[next].getGeneName() + ";after_ID=" + g[i].getGeneId() + ";after_gene_id=" + g[i].getGeneName() + "\n");
                                wr.write("chr" + g[next].getChr() + "\tPRED\tCON\t" + (enddown + 1) + "\t" + (startup - 1) + "\t.\t-\t.\tID=" + g[next].getGeneId() + ";gene_id=" + g[next].getGeneName() + ";after_ID=" + g[i].getGeneId() + ";after_gene_id=" + g[i].getGeneName() + "\n");
                                wr2.write(g[next].getGeneId() + "\t" + g[i].getGeneId() + "\t" + g[i].getChr() + "\t" + g[next].getGeneStart() + "\t" + g[i].getGeneEnd() + "\n");
                                wr3.write(g[i].getGeneId() + "\t" + g[next].getGeneId() + "\t" + g[i].getChr() + "\t" + strand + "\t" + g[i].getGeneStart() + "\t" + g[i].getGeneEnd() + "\t" + g[next].getGeneStart() + "\t" + g[next].getGeneEnd() + "\n");
                            } else
                            {

                                wr.write("chr" + g[i].getChr() + "\tPRED\tDOWN\t" + startdown + "\t" + enddown + "\t.\t+\t.\tID=" + g[i].getGeneId() + ";gene_id=" + g[i].getGeneName() + ";after_ID=" + g[next].getGeneId() + ";after_gene_id=" + g[next].getGeneName() + "\n");
                                wr.write("chr" + g[next].getChr() + "\tPRED\tUP\t" + startup + "\t" + endup + "\t.\t+\t.\tID=" + g[next].getGeneId() + ";gene_id=" + g[next].getGeneName() + ";before_ID=" + g[i].getGeneId() + ";before_gene_id=" + g[i].getGeneName() + "\n");
                                wr.write("chr" + g[i].getChr() + "\tPRED\tCON\t" + (enddown + 1) + "\t" + (startup - 1) + "\t.\t+\t.\tID=" + g[i].getGeneId() + ";gene_id=" + g[i].getGeneName() + ";after_ID=" + g[next].getGeneId() + ";after_gene_id=" + g[next].getGeneName() + "\n");
                                wr2.write(g[i].getGeneId() + "\t" + g[next].getGeneId() + "\t" + g[i].getChr() + "\t" + g[i].getGeneEnd() + "\t" + g[next].getGeneStart() + "\n");
                                wr3.write(g[i].getGeneId() + "\t" + g[next].getGeneId() + "\t" + g[i].getChr() + "\t" + strand + "\t" + g[i].getGeneStart() + "\t" + g[i].getGeneEnd() + "\t" + g[next].getGeneStart() + "\t" + g[next].getGeneEnd() + "\n");
                            }
                        }

                    }

                   
                    if (i == 0 || (!g[i].getChr().equals(g[i - 1].getChr())))
                    {

                        int endup = g[i].getGeneStart() - 1;
                        int startup = Math.max(1, g[i].getGeneStart() - maxLength);

                        if (strand < 0)
                        {
                            wr.write("chr" + g[i].getChr() + "\tPRED\tDOWN\t" + startup + "\t" + endup + "\t.\t-\t.\tID=" + g[i].getGeneId() + ";gene_id=" + g[i].getGeneName() + ";after_ID=NA" + "\n");

                        } else
                        {
                            wr.write("chr" + g[i].getChr() + "\tPRED\tUP\t" + startup + "\t" + endup + "\t.\t+\t.\tID=" + g[i].getGeneId() + ";gene_id=" + g[i].getGeneName() + ";before_ID=NA" + "\n");

                        }

                    }

                }

            }
            wr.close();
            wr2.close();
            wr3.close();

        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void getConnections(String outfile, int maxLength, int minLength)
    {
        try
        {
            GeneInformation[] g = this.genes.values().toArray(new GeneInformation[0]);
            Arrays.sort(g);

            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
            int rightMostEnd = -1;
            String chr = null;

            GeneInformation rightMostGene = null;
            GeneInformation gene;

            for (int i = 0; i < g.length; i++)
            {
                
                if (chr != null && (i + 1 < g.length && !g[i].getChr().equals(g[i + 1].getChr())) || (i + 1 == g.length))
                {
                    if (rightMostGene.getGeneEnd() >= g[i].getGeneEnd())
                    {
                        gene = rightMostGene;
                    } else
                    {
                        gene = g[i];
                    }
                    int enddown = gene.getGeneEnd() + maxLength;
                    int startdown = gene.getGeneEnd() + 1;

                    if (gene.getStrand() < 0)
                    {
                        wr.write("chr" + gene.getChr() + "\tPRED\tUP\t" + startdown + "\t" + enddown + "\t" + gene.getGeneId() + "\t" + gene.getGeneName() + "\tNA\tNA" + "\n");

                    } else
                    {
                        wr.write("chr" + gene.getChr() + "\tPRED\tDOWN\t" + startdown + "\t" + enddown + "\t" + gene.getGeneId() + "\t" + gene.getGeneName() + "\tNA\tNA" + "\n");

                    }
                }

                if (chr == null || !chr.equals(g[i].getChr()) || g[i].getGeneStart() > rightMostEnd)
                {

                    if (g[i].getGeneEnd() > rightMostEnd || !chr.equals(g[i].getChr()))
                    {
                        rightMostEnd = g[i].getGeneEnd();
                        rightMostGene = g[i];
                        chr = g[i].getChr();
                    }
                    wr.write(g[i].getUTRSegments());

                    int next = i + 1;

                    while (next < g.length && g[next].getGeneStart() <= g[i].getGeneEnd() && g[next].getGeneEnd() <= g[i].getGeneEnd() && g[next].getChr().equals(g[i].getChr()))
                    {
                        next++;

                    }

                    if (next < g.length && g[next].getChr().equals(g[i].getChr()))
                    {

                        if (g[next].getGeneStart() > g[i].getGeneEnd() + minLength)
                        {
                            int intervalLength = g[next].getGeneStart() - g[i].getGeneEnd() - 1;
                            int updownsize = Math.min(intervalLength / 3, maxLength);
                            int enddown = g[i].getGeneEnd() + updownsize;
                            int startdown = g[i].getGeneEnd() + 1;

                            int endup = g[next].getGeneStart() - 1;
                            int startup = g[next].getGeneStart() - updownsize;

                            if (g[i].getStrand() < 0)
                            {
                                wr.write("chr" + g[i].getChr() + "\tPRED\tUP\t" + startdown + "\t" + enddown + "\t" + g[i].getGeneId() + "\t" + g[i].getGeneName() + "\t" + g[next].getGeneId() + "\t" + g[next].getGeneName() + "\n");

                            } else
                            {
                                wr.write("chr" + g[i].getChr() + "\tPRED\tDOWN\t" + startdown + "\t" + enddown + "\t" + g[i].getGeneId() + "\t" + g[i].getGeneName() + "\t" + g[next].getGeneId() + "\t" + g[next].getGeneName() + "\n");
                            }

                            if (g[next].getStrand() < 0)
                            {
                                wr.write("chr" + g[next].getChr() + "\tPRED\tDOWN\t" + startup + "\t" + endup + "\t" + g[next].getGeneId() + "\t" + g[next].getGeneName() + "\t" + g[i].getGeneId() + "\t" + g[i].getGeneName() + "\n");

                            } else
                            {

                                wr.write("chr" + g[next].getChr() + "\tPRED\tUP\t" + startup + "\t" + endup + "\t" + g[next].getGeneId() + "\t" + g[next].getGeneName() + "\t" + g[i].getGeneId() + "\t" + g[i].getGeneName() + "\n");
                            }
                        }

                    }

                    if (i == 0 || (!g[i].getChr().equals(g[i - 1].getChr())))
                    {

                        int endup = g[i].getGeneStart() - 1;
                        int startup = Math.max(1, g[i].getGeneStart() - maxLength);

                        if (g[i].getStrand() < 0)
                        {
                            wr.write("chr" + g[i].getChr() + "\tPRED\tDOWN\t" + startup + "\t" + endup + "\t" + g[i].getGeneId() + "\t" + g[i].getGeneName() + "\tNA\tNA" + "\n");

                        } else
                        {
                            wr.write("chr" + g[i].getChr() + "\tPRED\tUP\t" + startup + "\t" + endup + "\t" + g[i].getGeneId() + "\t" + g[i].getGeneName() + "\tNA\tNA" + "\n");

                        }

                    }

                }

            }
            wr.close();

        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void getDistanceDownstream(String outfile, int maxDist)
    {
        try
        {
            GeneInformation[] g = this.genes.values().toArray(new GeneInformation[0]);
            Arrays.sort(g);

            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
            int rightMostEnd = -1;
            int leftMostEnd = -1;
            String chr = null;

            boolean done;

            int next;
            int dist;
            GeneInformation previous=null;
            wr.write("GeneID\tGeneSymbol\tStrand\tDistance\tDownstreamGeneID\tDownstreamGeneSymbol\tDownstreamGeneStrand\n");
            for (int i = 0; i < g.length; i++)
            {
                done = false;
                if (g[i].getStrand() > 0)
                {
                    if (chr != null && chr.equals(g[i].getChr()))
                    {
                        if (rightMostEnd > g[i].getGeneEnd())
                        {
                            dist=0;
                           // wr.write(g[i] + "\t" + dist + "\t" + previous.getGeneId() + "\t" + previous.getGeneName() + "\n");
                            wr.write(g[i].getGeneId()+"\t"+g[i].getGeneName()+"\t"+g[i].getStrand() + "\t" + dist + "\t" + previous.getGeneId() + "\t" + previous.getGeneName() + "\t" + previous.getStrand()+ "\n");
                            done = true;
                        }

                    }
                    if (!done)
                    {
                        next = i + 1;
                        while (next < g.length && g[next].getChr().equals(g[i].getChr()) && g[next].getGeneEnd() <= g[i].getGeneEnd())
                        {
                            next++;
                        }
                        if (next >= g.length || !g[next].getChr().equals(g[i].getChr()))
                        {
                            wr.write(g[i].getGeneId()+"\t"+g[i].getGeneName()+"\t"+g[i].getStrand() + "\t" + maxDist + "\tNA\tNA\tNA\n");
                        } else
                        {
                            if (g[next].getGeneStart() <= g[i].getGeneEnd())
                            {
                                wr.write(g[i].getGeneId()+"\t"+g[i].getGeneName()+"\t"+g[i].getStrand() + "\t0\t" + g[next].getGeneId() + "\t" + g[next].getGeneName() + "\t" + g[next].getStrand() + "\n");
                            } else
                            {
                                dist = g[next].getGeneStart() - g[i].getGeneEnd() + 1;
                                wr.write(g[i].getGeneId()+"\t"+g[i].getGeneName()+"\t"+g[i].getStrand() + "\t" + dist + "\t" + g[next].getGeneId() + "\t" + g[next].getGeneName() + "\t" + g[next].getStrand() + "\n");
                            }
                        }

                    }

                } else
                {
                    if (chr != null && chr.equals(g[i].getChr()))
                    {
                        if (rightMostEnd > g[i].getGeneStart())
                        {
                           
                            dist=0;
                        }
                        else{
                             dist = g[i].getGeneStart() - previous.getGeneEnd() + 1;
                            
                        }
                         wr.write(g[i].getGeneId()+"\t"+g[i].getGeneName()+"\t"+g[i].getStrand() + "\t" + dist + "\t" + previous.getGeneId() + "\t" + previous.getGeneName() + "\t" + previous.getStrand() + "\n");
                    }
                    else{
                         wr.write(g[i].getGeneId()+"\t"+g[i].getGeneName()+"\t"+g[i].getStrand() + "\t" + maxDist + "\tNA\tNA\tNA\n");
                    }

                }

                if (chr != null && chr.equals(g[i].getChr()))
                {
                    if (g[i].getGeneEnd() > rightMostEnd)
                    {
                        rightMostEnd = g[i].getGeneEnd();
                        previous=g[i];
                    }

                } else
                {
                    rightMostEnd = g[i].getGeneEnd();
                    previous=g[i];

                }

                chr = g[i].getChr();

            }
            wr.close();

        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void makeReadoutSAFFile(String outfile, int maxLength, int minLength, int strand, int windowlength, String chrSizesFile)
    {
        try
        {
            GeneInformation[] g = this.genes.values().toArray(new GeneInformation[0]);
            Arrays.sort(g);

            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
            int rightMostEnd = -1;
            String chr = null;

            HashMap<String, Integer> chrSizes = new HashMap<>();

            BufferedReader rd = new BufferedReader(new FileReader(chrSizesFile));
            String inline = rd.readLine();
            String[] content;

            while (inline != null)
            {
                content = inline.split("\t");
                chrSizes.put(content[0], Integer.parseInt(content[1]));

                inline = rd.readLine();
            }
            rd.close();

            String chrName;

            for (int i = 0; i < g.length; i++)
            {

                //if(i == 0 || g[i - 1].getGeneEnd() < g[i].getGeneEnd() || !g[i-1].getChr().equals(g[i].getChr()))
                if (chr == null || !chr.equals(g[i].getChr()) || g[i].getGeneStart() > rightMostEnd)
                {

                    if (g[i].getGeneEnd() > rightMostEnd || !chr.equals(g[i].getChr()))
                    {
                        rightMostEnd = g[i].getGeneEnd();
                        chr = g[i].getChr();
                    }

                    int next = i + 1;

                    while (next < g.length && g[next].getGeneStart() <= g[i].getGeneEnd() && g[next].getGeneEnd() <= g[i].getGeneEnd() && g[next].getChr().equals(g[i].getChr()))
                    {
                        next++;
                    }

                    chrName = g[i].getChr();
                    chrName = "chr" + chrName;
                    if (chrName.equals("chrMT"))
                    {
                        chrName = "chrM";
                    }

                    if (next < g.length && g[next].getChr().equals(g[i].getChr()) && g[next].getGeneStart() > g[i].getGeneEnd() + minLength)
                    {

                        if (strand > 0)
                        {
                            int intervalStart = g[i].getGeneEnd() + 1;
                            int intervalEnd = Math.min(g[next].getGeneStart(), intervalStart + maxLength) - 1;

                            int start = intervalStart;
                            while (start + windowlength - 1 <= intervalEnd)
                            {

                                wr.write(g[i].getGeneId() + "\t" + chrName + "\t" + start + "\t" + (start + windowlength - 1) + "\t+\n");
                                start += windowlength;
                            }

                            if (start <= intervalEnd)
                            {
                                wr.write(g[i].getGeneId() + "\t" + chrName + "\t" + start + "\t" + intervalEnd + "\t+\n");
                            }

                        } else
                        {

                            int intervalEnd = g[next].getGeneStart() - 1;
                            int intervalStart = Math.max(g[i].getGeneEnd(), intervalEnd - maxLength) + 1;

                            int start = intervalEnd;
                            while (start - windowlength + 1 >= intervalStart)
                            {

                                wr.write(g[next].getGeneId() + "\t" + chrName + "\t" + (start - windowlength + 1) + "\t" + (start) + "\t-\n");
                                start -= windowlength;
                            }

                            if (start >= intervalStart)
                            {
                                wr.write(g[next].getGeneId() + "\t" + chrName + "\t" + intervalStart + "\t" + start + "\t-\n");
                            }

                        }

                    }

                    
                    if (strand > 0 & (next >= g.length || !g[i].getChr().equals(g[next].getChr())))
                    {

                        int chrLength = chrSizes.get(chrName);

                        int intervalStart = g[i].getGeneEnd() + 1;
                        int intervalEnd = Math.min(intervalStart + maxLength - 1, chrLength);

                        int start = intervalStart;
                        while (start + windowlength - 1 <= intervalEnd)
                        {

                            wr.write(g[i].getGeneId() + "\t" + chrName + "\t" + start + "\t" + (start + windowlength - 1) + "\t+\n");
                            start += windowlength;
                        }

                        if (start <= intervalEnd)
                        {
                            wr.write(g[i].getGeneId() + "\t" + chrName + "\t" + start + "\t" + intervalEnd + "\t+\n");
                        }

                    }
                    if (strand < 0 && (i == 0 || (next < g.length && !g[i].getChr().equals(g[next].getChr()))))
                    {
                        if (i == 0)
                        {
                            next = i;
                        }

                        chrName = g[next].getChr();
                        chrName = "chr" + chrName;
                        if (chrName.equals("chrMT"))
                        {
                            chrName = "chrM";
                        }
                        int intervalEnd = g[next].getGeneStart() - 1;
                        int intervalStart = Math.max(0, intervalEnd - maxLength) + 1;

                        int start = intervalEnd;
                        while (start - windowlength + 1 >= intervalStart)
                        {

                            wr.write(g[next].getGeneId() + "\t" + chrName + "\t" + (start - windowlength + 1) + "\t" + (start) + "\t-\n");
                            start -= windowlength;
                        }

                        if (start >= intervalStart)
                        {
                            wr.write(g[next].getGeneId() + "\t" + chrName + "\t" + intervalStart + "\t" + start + "\t-\n");
                        }

                    }

                }

            }
            wr.close();

        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void makeReadinSAFFile(String outfile, int maxLength, int minLength, int strand, int windowlength, String chrSizesFile)
    {
        try
        {
            GeneInformation[] g = this.genes.values().toArray(new GeneInformation[0]);
            Arrays.sort(g);

            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
            int rightMostEnd = -1;
            String chr = null;

            HashMap<String, Integer> chrSizes = new HashMap<String, Integer>();

            BufferedReader rd = new BufferedReader(new FileReader(chrSizesFile));
            String inline = rd.readLine();
            String[] content;

            while (inline != null)
            {
                content = inline.split("\t");
                chrSizes.put(content[0], Integer.parseInt(content[1]));

                inline = rd.readLine();
            }
            rd.close();

            String chrName;

            for (int i = 0; i < g.length; i++)
            {

                
                if (chr == null || !chr.equals(g[i].getChr()) || g[i].getGeneStart() > rightMostEnd)
                {

                    if (g[i].getGeneEnd() > rightMostEnd || !chr.equals(g[i].getChr()))
                    {
                        rightMostEnd = g[i].getGeneEnd();
                        chr = g[i].getChr();
                    }

                    int next = i + 1;

                    while (next < g.length && g[next].getGeneStart() <= g[i].getGeneEnd() && g[next].getGeneEnd() <= g[i].getGeneEnd() && g[next].getChr().equals(g[i].getChr()))
                    {
                        next++;
                    }

                    chrName = g[i].getChr();
                    chrName = "chr" + chrName;
                    if (chrName.equals("chrMT"))
                    {
                        chrName = "chrM";
                    }

                    if (next < g.length && g[next].getChr().equals(g[i].getChr()) && g[next].getGeneStart() > g[i].getGeneEnd() + minLength)
                    {

                        if (strand > 0)
                        {
                            int intervalStart = g[i].getGeneEnd() + 1;
                            int intervalEnd = Math.min(g[next].getGeneStart(), intervalStart + maxLength) - 1;

                            int start = intervalStart;
                            while (start + windowlength - 1 <= intervalEnd)
                            {

                                wr.write(g[i].getGeneId() + "\t" + chrName + "\t" + start + "\t" + (start + windowlength - 1) + "\t+\n");
                                start += windowlength;
                            }

                            if (start <= intervalEnd)
                            {
                                wr.write(g[i].getGeneId() + "\t" + chrName + "\t" + start + "\t" + intervalEnd + "\t+\n");
                            }

                        } else
                        {

                            int intervalEnd = g[next].getGeneStart() - 1;
                            int intervalStart = Math.max(g[i].getGeneEnd(), intervalEnd - maxLength) + 1;

                            int start = intervalEnd;
                            while (start - windowlength + 1 >= intervalStart)
                            {

                                wr.write(g[next].getGeneId() + "\t" + chrName + "\t" + (start - windowlength + 1) + "\t" + (start) + "\t-\n");
                                start -= windowlength;
                            }

                            if (start >= intervalStart)
                            {
                                wr.write(g[next].getGeneId() + "\t" + chrName + "\t" + intervalStart + "\t" + start + "\t-\n");
                            }

                        }

                    }

                   
                    if (strand > 0 & (next >= g.length || !g[i].getChr().equals(g[next].getChr())))
                    {

                        int chrLength = chrSizes.get(chrName);

                        int intervalStart = g[i].getGeneEnd() + 1;
                        int intervalEnd = Math.min(intervalStart + maxLength - 1, chrLength);

                        int start = intervalStart;
                        while (start + windowlength - 1 <= intervalEnd)
                        {

                            wr.write(g[i].getGeneId() + "\t" + chrName + "\t" + start + "\t" + (start + windowlength - 1) + "\t+\n");
                            start += windowlength;
                        }

                        if (start <= intervalEnd)
                        {
                            wr.write(g[i].getGeneId() + "\t" + chrName + "\t" + start + "\t" + intervalEnd + "\t+\n");
                        }

                    }
                    if (strand < 0 && (i == 0 || (next < g.length && !g[i].getChr().equals(g[next].getChr()))))
                    {
                        if (i == 0)
                        {
                            next = i;
                        }

                        chrName = g[next].getChr();
                        chrName = "chr" + chrName;
                        if (chrName.equals("chrMT"))
                        {
                            chrName = "chrM";
                        }
                        int intervalEnd = g[next].getGeneStart() - 1;
                        int intervalStart = Math.max(0, intervalEnd - maxLength) + 1;

                        int start = intervalEnd;
                        while (start - windowlength + 1 >= intervalStart)
                        {

                            wr.write(g[next].getGeneId() + "\t" + chrName + "\t" + (start - windowlength + 1) + "\t" + (start) + "\t-\n");
                            start -= windowlength;
                        }

                        if (start >= intervalStart)
                        {
                            wr.write(g[next].getGeneId() + "\t" + chrName + "\t" + intervalStart + "\t" + start + "\t-\n");
                        }

                    }

                }

            }
            wr.close();

        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void countJunctionReads(String infile, String outfile, int strand)
    {
        GeneInformation[] g = new GeneInformation[this.genes.size()];
        int i = 0;
        for (GeneInformation gene : genes.values())
        {
            if (gene.getStrand() * strand > 0)
            {
                g[i] = gene;
                i++;
            }
        }
        GeneInformation[] tmp = new GeneInformation[i];
        System.arraycopy(g, 0, tmp, 0, i);
        g = tmp;
        Arrays.sort(g);
        int pos = 0;

        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(infile));
            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));

            String inline = rd.readLine();
            String[] content;
            String chr;
            String cigar;
            Pattern pattern = Pattern.compile("(\\d+)M(\\d+)N(\\d+)M");
            Matcher m;
            int length;
            int start;
            int end;
            int maxlength;
            int bestGene;
            int[] counts = new int[g.length];
            while (inline != null)
            {
                content = inline.split("\t");
                cigar = content[5];
                m = pattern.matcher(cigar);
                if (m.matches())
                {
                    length = Integer.parseInt(m.group(1)) + Integer.parseInt(m.group(2)) + Integer.parseInt(m.group(3));

                    chr = content[2];
                    if (chr.startsWith("chr"))
                    {
                        chr = chr.substring(3);
                    }
                    while (pos < g.length && !g[pos].getChr().equals(chr))
                    {
                        pos++;
                    }
                    start = Integer.parseInt(content[3]);

                    end = start + length - 1;

                    while (pos < g.length && g[pos].getGeneEnd() < start)
                    {
                        pos++;
                    }
                    i = pos;
                    maxlength = 0;
                    bestGene = -1;
                    while (i < g.length && g[i].getGeneStart() <= start)
                    {
                        if (g[i].getGeneEnd() >= end)
                        {
                            length = g[i].getGeneEnd() - g[i].getGeneStart() + 1;
                            if (length > maxlength)
                            {
                                maxlength = length;
                                bestGene = i;

                            }
                        }
                        i++;
                    }
                    if (bestGene != -1)
                    {
                        counts[bestGene]++;

                    }

                }
                inline = rd.readLine();
            }
            rd.close();

            for (int j = 0; j < counts.length; j++)
            {
                wr.write(g[j].getGeneId() + "\t" + counts[j] + "\tchr" + g[j].getChr() + "\t" + g[j].getGeneStart() + "\t" + g[j].getGeneEnd() + "\t" + g[j].getStrand() + "\n");
            }

            wr.close();
        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void getRelevantRegions(int distToNextGene, int downstreamLength, int minLength, String outfile)
    {

        try
        {
            GeneInformation[] g = genes.values().toArray(new GeneInformation[0]);
            int i = 0;
            for (GeneInformation gene : genes.values())
            {
                g[i] = gene;
                i++;
            }
            Arrays.sort(g);

            BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));

            GeneInformation before = null;
            GeneInformation after = null;
            int endpos;
            boolean containedWithinOtherGene;
            for (int j = 0; j < g.length; j++)
            {
                if (g[j].getChr().length() <= 2)
                {
                    
                    {

                        if (g[j].getStrand() < 0)
                        {

                            containedWithinOtherGene = false;
                            if (j != 0 && g[j - 1].getChr().equals(g[j].getChr()))
                            {

                                before = g[j - 1];
                                int k = j - 2;

                                if (before.getGeneEnd() >= g[j].getGeneEnd() && before.getStrand() == g[j].getStrand())
                                {
                                    containedWithinOtherGene = true;

                                    // System.out.println(g[j]+"\n"+before+"\n");
                                }
                                while (k > 0 && j - k <= 10 && !containedWithinOtherGene)
                                {

                                    if (g[k].getGeneEnd() > before.getGeneEnd() || before.getStrand() != g[j].getStrand())
                                    {
                                        if (g[k].getStrand() == g[j].getStrand())
                                        {
                                            before = g[k];

                                            if (before.getGeneEnd() >= g[j].getGeneEnd() && before.getStrand() == g[j].getStrand())
                                            {
                                                containedWithinOtherGene = true;
                                                //System.out.println(g[j]+"\n"+before+"\n");
                                            }
                                        }
                                    }
                                    k--;
                                }
                                endpos = before.getGeneEnd() + distToNextGene;
                                endpos = Math.max(endpos, g[j].getGeneStart() - downstreamLength);
                                if (!containedWithinOtherGene && before.getGeneEnd() < g[j].getGeneStart() - 1)
                                {
                                    wr.write("chr" + g[j].getChr() + "\tPRED\tCON\t" + (before.getGeneEnd() + 1) + "\t" + (g[j].getGeneStart() - 1) + "\t.\t-\t.\tID=" + before.getGeneId() + "-" + g[j].getGeneId() + ";gene_id=" + before.getGeneName() + "-" + g[j].getGeneName() + ";note=CON\n");

                                }

                            } else
                            {
                                endpos = g[j].getGeneStart() - downstreamLength;
                            }
                            endpos = Math.max(1, endpos);

                            if (!containedWithinOtherGene)
                            {
                                if (g[j].hasUTR())
                                {
                                    wr.write("chr" + g[j].getChr() + "\tPRED\tUTR\t" + g[j].getGeneStart() + "\t" + g[j].getLast3PrimeUTREnd() + "\t.\t-\t.\tID=" + g[j].getGeneId() + ";gene_id=" + g[j].getGeneName() + ";note=3UTR\n");
                                }
                                if (g[j].getGeneStart() - endpos >= minLength)
                                {
                                    wr.write("chr" + g[j].getChr() + "\tPRED\tDOWN\t" + endpos + "\t" + (g[j].getGeneStart() - 1) + "\t.\t-\t.\tID=" + g[j].getGeneId() + ";gene_id=" + g[j].getGeneName() + ";note=DOWN\n");

                                }

                                if (j != g.length - 1 && g[j + 1].getChr().equals(g[j].getChr()))
                                {

                                    int k = j + 1;
                                    while (k < g.length && (g[k].getGeneEnd() <= g[j].getGeneEnd() || g[k].getStrand() != g[j].getStrand()))
                                    {
                                        k++;
                                    }

                                    if (k < g.length && g[k].getChr().equals(g[j].getChr()))
                                    {
                                        after = g[k];
                                        endpos = after.getGeneStart() - distToNextGene;
                                        endpos = Math.min(endpos, g[j].getGeneEnd() + downstreamLength);

                                    } else
                                    {
                                        endpos = g[j].getGeneEnd() + downstreamLength;
                                    }
                                } else
                                {
                                    endpos = g[j].getGeneEnd() + downstreamLength;
                                }

                                if (endpos - g[j].getGeneEnd() >= minLength)
                                {

                                    wr.write("chr" + g[j].getChr() + "\tPRED\tUP\t" + (g[j].getGeneEnd() + 1) + "\t" + endpos + "\t.\t-\t.\tID=" + g[j].getGeneId() + ";gene_id=" + g[j].getGeneName() + ";note=UP\n");
                                }

                            }

                        } else
                        {
                            containedWithinOtherGene = false;

                            if (j != 0 && g[j - 1].getChr().equals(g[j].getChr()))
                            {
                                before = g[j - 1];
                                int k = j - 2;
                                if (before.getGeneEnd() >= g[j].getGeneEnd() && before.getStrand() == g[j].getStrand())
                                {
                                    containedWithinOtherGene = true;
                                    

                                }

                                while (k > 0 && j - k <= 10)
                                {
                                    if (g[k].getGeneEnd() > before.getGeneEnd() || before.getStrand() != g[j].getStrand())
                                    {
                                        if (g[k].getStrand() == g[j].getStrand())
                                        {
                                            before = g[k];
                                            if (before.getGeneEnd() >= g[j].getGeneEnd() && before.getStrand() == g[j].getStrand())
                                            {
                                                containedWithinOtherGene = true;
                                               

                                            }
                                        }
                                    }
                                    k--;
                                }
                                endpos = before.getGeneEnd() + distToNextGene;
                                endpos = Math.max(endpos, g[j].getGeneStart() - downstreamLength);

                            } else
                            {
                                endpos = g[j].getGeneStart() - downstreamLength;

                            }
                            endpos = Math.max(1, endpos);
                            if (!containedWithinOtherGene)
                            {
                                if (g[j].getGeneStart() - endpos >= minLength)
                                {
                                    wr.write("chr" + g[j].getChr() + "\tPRED\tUP\t" + endpos + "\t" + (g[j].getGeneStart() - 1) + "\t.\t+\t.\tID=" + g[j].getGeneId() + ";gene_id=" + g[j].getGeneName() + ";note=UP\n");

                                }

                                if (j != g.length - 1 && g[j + 1].getChr().equals(g[j].getChr()))
                                {

                                    int k = j + 1;
                                    while (k < g.length && (g[k].getGeneEnd() <= g[j].getGeneEnd() || g[k].getStrand() != g[j].getStrand()))
                                    {
                                        k++;
                                    }

                                    if (k < g.length && g[k].getChr().equals(g[j].getChr()))
                                    {
                                        after = g[k];

                                        endpos = after.getGeneStart() - distToNextGene;
                                        endpos = Math.min(endpos, g[j].getGeneEnd() + downstreamLength);
                                        if (g[j].getGeneEnd() < after.getGeneStart() - 1)
                                        {
                                            wr.write("chr" + g[j].getChr() + "\tPRED\tCON\t" + (g[j].getGeneEnd() + 1) + "\t" + (after.getGeneStart() - 1) + "\t.\t+\t.\tID=" + g[j].getGeneId() + "-" + after.getGeneId() + ";gene_id=" + g[j].getGeneName() + "-" + after.getGeneName() + ";note=CON\n");
                                        }

                                    } else
                                    {
                                        endpos = g[j].getGeneEnd() + downstreamLength;
                                    }
                                } else
                                {
                                    endpos = g[j].getGeneEnd() + downstreamLength;
                                }

                                if (g[j].hasUTR())
                                {
                                    wr.write("chr" + g[j].getChr() + "\tPRED\tUTR\t" + g[j].getLast3PrimeUTRStart() + "\t" + g[j].getGeneEnd() + "\t.\t+\t.\tID=" + g[j].getGeneId() + ";gene_id=" + g[j].getGeneName() + ";note=3UTR\n");
                                }
                                if (endpos - g[j].getGeneEnd() >= minLength)
                                {

                                    wr.write("chr" + g[j].getChr() + "\tPRED\tDOWN\t" + (g[j].getGeneEnd() + 1) + "\t" + endpos + "\t.\t+\t.\tID=" + g[j].getGeneId() + ";gene_id=" + g[j].getGeneName() + ";note=DOWN\n");
                                }

                            }

                        }

                    }
                }
            }

            wr.close();
        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void printGeneInformation()
    {

        ArrayList<GeneInformation> g = new ArrayList<GeneInformation>();
        for (GeneInformation gene : genes.values())
        {
            g.add(gene);
        }
        Collections.sort(g);

        for (GeneInformation gene : g)
        {
            System.out.println(gene);
        }

    }

    public void analyzeConnections(String infile, String clusterFile, int minClusterUp, int minClusterDown)
    {
        try
        {
            BufferedReader rd = new BufferedReader(new FileReader(clusterFile));
            String inline = rd.readLine();
            String[] content;
            HashMap<String, Integer> clustersUp = new HashMap<>();
            HashMap<String, Integer> clustersDown = new HashMap<>();
            int clusterUp;
            int clusterDown;

            while (inline != null)
            {
                content = inline.split("\t");
                if (!content[1].equals("NA"))
                {
                    clusterUp = Integer.parseInt(content[1]);
                    clustersUp.put(content[0], clusterUp);
                }
                if (!content[2].equals("NA"))
                {
                    clusterDown = Integer.parseInt(content[2]);
                    clustersDown.put(content[0], clusterDown);
                }

                inline = rd.readLine();
            }
            rd.close();

            rd = new BufferedReader(new FileReader(infile));
            inline = rd.readLine();

            Pattern pattern = Pattern.compile("ID=(.*);gene_id=(.*);after_ID=(.*);after_gene_id=(.*)");

            HashMap<String, String> connections = new HashMap<>();
            Matcher m;
            String gene1;
            String gene2;
            String geneSymbol1;
            String geneSymbol2;
            HashSet<String> starts = new HashSet<>();
            while (inline != null)
            {
                content = inline.split("\t");
                if (content[2].equals("CON"))
                {
                    m = pattern.matcher(inline);
                    if (m.find())
                    {
                        gene1 = m.group(1);
                        gene2 = m.group(3);
                        geneSymbol1 = m.group(2);
                        geneSymbol2 = m.group(4);

                        if (clustersDown.containsKey(gene1) && clustersDown.get(gene1) >= minClusterDown)
                        {
                            if (clustersUp.containsKey(gene2))
                            {
                                if (clustersUp.get(gene2) >= minClusterUp)
                                {
                                    connections.put(geneSymbol1, geneSymbol2);
                                } else
                                {
                                    connections.put(geneSymbol1, "STOP NO INREAD");
                                }
                            } else
                            {
                                connections.put(geneSymbol1, "STOP NO EXP");
                            }
                           

                            if (clustersUp.containsKey(gene1) && clustersUp.get(gene1) < minClusterUp)
                            {
                                starts.add(geneSymbol1);

                            }
                        }

                    }
                }
                inline = rd.readLine();
            }
            rd.close();

           
            String s;
            for (String k : starts)
            {
                s = k;

                while (s != null)
                {
                    System.out.print(s + "\t");
                    s = connections.get(s);

                }
                System.out.println();

            }
        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public GeneInformation[] loadPTT(String infile, int minLength)
    {
        ArrayList<GeneInformation> ptts = new ArrayList<>();
        try
        {

            BufferedReader rd = new BufferedReader(new FileReader(infile));
            String inline = rd.readLine();
            String[] content;
            String chr;
            int start;
            int end;

            GeneInformation ptt;
            while (inline != null)
            {
                content = inline.split("\t");
                chr = content[0].substring(3);
                start = Integer.parseInt(content[1]);
                end = Integer.parseInt(content[2]);
                if (end - start + 1 >= minLength)
                {
                    ptt = new GeneInformation();
                    ptt.setGeneName("PTT");
                    ptt.setChr(chr);
                    ptt.setGeneStart(start);
                    ptt.setGeneEnd(end);

                    ptts.add(ptt);
                }

                inline = rd.readLine();
            }
            rd.close();

        } catch (IOException ex)
        {
            Logger.getLogger(UTRCalculator.class.getName()).log(Level.SEVERE, null, ex);
        }
        GeneInformation[] ptts_sorted = ptts.toArray(new GeneInformation[0]);
        Arrays.sort(ptts_sorted);
        return ptts_sorted;
    }

    
}
