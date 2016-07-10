package nl.peterbloem.motive.exec;

import static nl.peterbloem.kit.Functions.log2;
import static nl.peterbloem.kit.Series.series;
import static org.apache.commons.math3.util.ArithmeticUtils.binomialCoefficientLog;
import static org.nodes.models.USequenceEstimator.CIMethod;
import static org.nodes.models.USequenceEstimator.CIType;
import static org.nodes.motifs.MotifCompressor.exDegree;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.json.JSONObject;
import org.nodes.DGraph;
import org.nodes.DLink;
import org.nodes.DNode;
import org.nodes.Graph;
import org.nodes.Graphs;
import org.nodes.Subgraph;
import org.nodes.TGraph;
import org.nodes.TLink;
import org.nodes.UGraph;
import org.nodes.ULink;
import org.nodes.UNode;
import org.nodes.algorithms.Nauty;
import org.nodes.compression.BinomialCompressor;
import org.nodes.compression.EdgeListCompressor;
import org.nodes.compression.NeighborListCompressor;
import org.nodes.data.Data;
import org.nodes.models.DSequenceEstimator;
import org.nodes.models.DegreeSequenceModel.Prior;
import org.nodes.models.ERSimpleModel;
import org.nodes.models.EdgeListModel;
import org.nodes.models.USequenceEstimator;
import org.nodes.motifs.MotifCompressor;
import org.nodes.random.RandomGraphs;
import org.nodes.random.SimpleSubgraphGenerator;
import org.nodes.util.bootstrap.BCaCI;
import org.nodes.util.bootstrap.LogBCaCI;
import org.nodes.util.bootstrap.LogNormalCI;

import nl.peterbloem.kit.FileIO;
import nl.peterbloem.kit.Functions;
import nl.peterbloem.kit.Generator;
import nl.peterbloem.kit.Generators;
import nl.peterbloem.kit.Global;
import nl.peterbloem.kit.OnlineModel;
import nl.peterbloem.kit.Order;
import nl.peterbloem.kit.Series;
import nl.peterbloem.motive.DPlainMotifExtractor;
import nl.peterbloem.motive.MotifModel;
import nl.peterbloem.motive.MotifSearchModel;
import nl.peterbloem.motive.UPlainMotifExtractor;

/**
 * Compares the code length under the motifs to that under a given null-model
 * 
 * For undirected data.
 * 
 * @author Peter
 */

public class CompareLarge
{	
	private static final int BS_SAMPLES = 10000;

	/** 
	 * Number of samples to take to find potential motifs 
	 */
	public int motifSamples = 1000000;
	/**
	 * Minimum motif size (inclusive)
	 */
	public int motifMinSize = 3;
	/**
	 * Maximum motif size (inclusive)
	 */
	public int motifMaxSize = 6;
	
	/**
	 * The dataset.
	 */
	public DGraph<String> data = null;
	/**
	 * Name of the dataset (to appear in the plot) 
	 */
	public String dataName = "<no dataset name given>";
	
	/**
	 * The maximum number of motifs to test
	 */
	public int maxMotifs = 100;
	
	/** 
	 * Minimum frequency for a motif to be considered
	 */
	public int minFreq = 2;

	/**
	 * Depth to which to search for which instances to consider (more is better, but slower, -1 is maximal depth always).
	 */
	public int searchDepth = -1;
	
	/**
	 * Whether to reset the DM model at every motif instance.
	 */
	private boolean resets = true;
	
	public void main() throws IOException
	{		
		nl.peterbloem.kit.Global.secureRandom(42);
		
		Global.log().info("Computing motif code lengths");
		

		DPlainMotifExtractor<String> ex 
		= new DPlainMotifExtractor<String>(
				(DGraph<String>)data, motifSamples, motifMinSize, motifMaxSize, minFreq);
	
		List<? extends DGraph<String>> subs = 
				new ArrayList<DGraph<String>>(ex.subgraphs());
		List<Double> frequencies = new ArrayList<Double>(subs.size());
		for(Graph<String> sub : subs)
			frequencies.add(ex.frequency((DGraph<String>)sub));
		
		List<List<List<Integer>>> occurrences = 
				new ArrayList<List<List<Integer>>>(subs.size());
		for(Graph<String> sub : subs)
			occurrences.add(ex.occurrences((DGraph<String>)sub));
	
		if(subs.size() > maxMotifs)
		{
			subs = new ArrayList<DGraph<String>>(subs.subList(0, maxMotifs));
			frequencies = new ArrayList<Double>(frequencies.subList(0, maxMotifs));
		}
					
		List<Double> factorsER = new ArrayList<Double>(subs.size());
		List<Double> factorsEL = new ArrayList<Double>(subs.size());
		List<Double> maxFactors = new ArrayList<Double>(subs.size());

		double baselineER = (new ERSimpleModel(false)).codelength(data);
		double baselineEL = (new EdgeListModel(Prior.ML)).codelength(data);
		
		for(int i : series(subs.size()))
		{
			DGraph<String> sub = subs.get(i);
			List<List<Integer>> occs = occurrences.get(i);
			
			Global.log().info("Analysing sub ("+ (i+1) +" of " + subs.size() + "): " + sub);
			Global.log().info("freq: " + frequencies.get(i));
			
			double max = Double.NEGATIVE_INFINITY;

			Global.log().info("null model: ER");
			{
				double sizeER = MotifSearchModel.sizeER(data, sub, occs, resets, searchDepth); 
				double factorER = baselineER - sizeER;
				factorsER.add(factorER);
				 
				Global.log().info("ER baseline: " + baselineER);
				Global.log().info("ER motif code: " + sizeER);
				Global.log().info("ER factor: " + factorER);
				
				max = Math.max(max, factorER);
			}

			Global.log().info("null model: EL");
			{
				double sizeEL = MotifSearchModel.sizeEL(data, sub, occs, resets, searchDepth); 
				double factorEL = baselineEL - sizeEL;
				factorsEL.add(factorEL);
			 
				Global.log().info("EL baseline: " + baselineEL);
				Global.log().info("EL motif code: " + sizeEL);
				Global.log().info("EL factor: " + factorEL);
				
				max = Math.max(max, factorEL);
			}

			maxFactors.add(max);
		}
		
		Comparator<Double> comp = Functions.natural();
		Functions.sort(
				factorsEL, Collections.reverseOrder(comp), 
				(List) frequencies,
				(List) factorsER, 
				(List) factorsEL, 
				(List) subs);
		
		File numbersFile = new File("numbers.csv");
		
		BufferedWriter numbersWriter = new BufferedWriter(new FileWriter(numbersFile));
		for(int i : series(subs.size()))
			numbersWriter.write(frequencies.get(i) + ", " + factorsER.get(i) + ", " + factorsEL.get(i) + "\n");		
		numbersWriter.close();

		int i = 0;
		for(Graph<String> sub : subs)
		{
			File graphFile = new File(String.format("motif.%03d.edgelist", i));
			Data.writeEdgeList(sub, graphFile);
			
			i++;
		}

		JSONObject obj = new JSONObject();
		obj.put("data", dataName);
		obj.put("directed", true);
		obj.put("baseline er", baselineER);
		obj.put("baseline el", baselineEL);
		Functions.write(obj.toString(), new File("metadata.json"));
		
		try
		{
			FileIO.python(new File("."), "scripts/plot.large.py");
		} catch (Exception e)
		{
			System.out.println("Failed to run plot script. " + e);
		}
	}
}
