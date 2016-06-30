package nl.peterbloem.motive.exec;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static nl.peterbloem.kit.Functions.log2Choose;
import static nl.peterbloem.kit.Functions.log2Factorial;
import static nl.peterbloem.kit.Functions.prefix;
import static nl.peterbloem.kit.OnlineModel.storeSequence;
import static nl.peterbloem.kit.OnlineModel.storeSequenceML;
import static nl.peterbloem.kit.Series.series;
import static org.nodes.Graphs.degrees;
import static org.nodes.models.USequenceEstimator.CIMethod;
import static org.nodes.models.USequenceEstimator.CIType;
import static org.nodes.motifs.MotifCompressor.MOTIF_SYMBOL;
import static org.nodes.motifs.MotifCompressor.exDegree;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.json.JSONObject;
import org.nodes.DGraph;
import org.nodes.DNode;
import org.nodes.Graph;
import org.nodes.Graphs;
import org.nodes.Link;
import org.nodes.Node;
import org.nodes.Subgraph;
import org.nodes.UGraph;
import org.nodes.ULink;
import org.nodes.UNode;
import org.nodes.algorithms.Nauty;
import org.nodes.clustering.ConnectionClusterer;
import org.nodes.compression.BinomialCompressor;
import org.nodes.compression.EdgeListCompressor;
import org.nodes.compression.NeighborListCompressor;
import org.nodes.data.Data;
import org.nodes.models.DSequenceEstimator;
import org.nodes.models.DegreeSequenceModel;
import org.nodes.models.DegreeSequenceModel.Margin;
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
import nl.peterbloem.kit.FrequencyModel;
import nl.peterbloem.kit.Functions;
import nl.peterbloem.kit.Generator;
import nl.peterbloem.kit.Generators;
import nl.peterbloem.kit.Global;
import nl.peterbloem.kit.OnlineModel;
import nl.peterbloem.kit.Order;
import nl.peterbloem.kit.Pair;
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
public class Compare
{
	/** 
	 * Number of samples to take to find potential motifs 
	 */
	public int motifSamples;
	
	/**
	 * Minimum motif size (inclusive)
	 */
	public int motifMinSize = 3;
	/**
	 * Maximum motif size (inclusive)
	 */
	public int motifMaxSize = 6;
	
	/**
	 * The number of samples to take for the DS model
	 */
	public int betaIterations = 50;
	
	/**
	 * The alpha to use in construction significance intervals for the DS model.
	 */
	public double betaAlpha = 0.05;

	/**
	 * The dataset.
	 */
	public Graph<String> data = null;
	
	/**
	 * Name of the dataset (to appear in the plot) 
	 */
	public String dataName = "";
	
	/**
	 * The maximum number of motifs to test
	 */
	public int maxMotifs = 100;
	
	/** 
	 * Minimum frequency for a motif to be considered
	 */
	public int minFreq = 2;
	
	/**
	 * Whether to transform the graph to a simple graph before motif extraction
	 */
	public boolean simplify = true;
	
	/**
	 * The depth to which to search when using the DS model
	 */
	public int betaSearchDepth = -1;
	
	/**
	 * Number of threads to use when sampling for the DS model.
	 */
	public static final int NUM_THREADS = Runtime.getRuntime().availableProcessors();;
	
	public static enum NullModel{ER, EDGELIST, BETA}
	
	boolean directed;

	/**
	 * Whether to reset the DM moel for each motif instance
	 */
	private boolean resets = true;
	
	public void main() throws IOException
	{
		Global.secureRandom(42);
		Global.log().info("Threads available: " +  NUM_THREADS);
		
		directed = data instanceof DGraph<?>;
		Global.log().info("Is data directed? : " + directed + " ("+data.getClass()+")");
		if(simplify)
		{
			if(directed)
				data = Graphs.toSimpleDGraph((DGraph<String>)data);
			else
				data = Graphs.toSimpleUGraph(data);
		}
		
		Global.log().info("data nodes: " + data.size());
		Global.log().info("data links: " + data.numLinks());
		
		List<Integer> degrees = Graphs.degrees(data);
		Collections.sort(degrees, Collections.reverseOrder());
		
		Global.log().info("Computing motif code lengths");
		
		List<? extends Graph<String>> subs;
		List<Double> frequencies;
		List<List<List<Integer>>> occurrences;

		if(directed)
		{
			DPlainMotifExtractor<String> ex 
			= new DPlainMotifExtractor<String>(
					(DGraph<String>)data, motifSamples, motifMinSize, motifMaxSize, minFreq);
		
			subs = new ArrayList<Graph<String>>(ex.subgraphs());
			
			frequencies = new ArrayList<Double>(subs.size());
			for(Graph<String> sub : subs)
				frequencies.add((double)ex.occurrences((DGraph<String>)sub).size());
			
			occurrences = new ArrayList<List<List<Integer>>>(subs.size());
			for(Graph<String> sub : subs)
				occurrences.add(ex.occurrences((DGraph<String>)sub));
		} else
		{	
			UPlainMotifExtractor<String> ex 
				= new UPlainMotifExtractor<String>(
						(UGraph<String>)data, motifSamples, motifMinSize, motifMaxSize, minFreq);
			
			subs = new ArrayList<Graph<String>>(ex.subgraphs());
			frequencies = new ArrayList<Double>(subs.size());
			for(Graph<String> sub : subs)
				frequencies.add((double)ex.occurrences((UGraph<String>)sub).size());
			
			occurrences = new ArrayList<List<List<Integer>>>(subs.size());
			for(Graph<String> sub : subs)
				occurrences.add(ex.occurrences((UGraph<String>)sub));
		}
		
		if(subs.size() > maxMotifs)
		{
			subs = new ArrayList<Graph<String>>(subs.subList(0, maxMotifs));
			frequencies = new ArrayList<Double>(frequencies.subList(0, maxMotifs));
		}
			
		List<Double> factorsER   = new ArrayList<Double>(subs.size());
		List<Double> factorsEL   = new ArrayList<Double>(subs.size());
		List<Double> factorsBeta = new ArrayList<Double>(subs.size());
		List<Double> maxFactors   =  new ArrayList<Double>(subs.size());
				
		double baselineER = new ERSimpleModel(true).codelength(data);
		double baselineEL = new EdgeListModel(Prior.ML).codelength(data);
		double baselineBeta = new DegreeSequenceModel(betaIterations, betaAlpha, Prior.ML, Margin.LOWERBOUND).codelength(data);
				
		for(int i : series(subs.size()))
		{
			Graph<String> sub = subs.get(i);
			List<List<Integer>> occs = occurrences.get(i);
			
			Global.log().info("Analysing sub ("+ (i+1) +" of " + subs.size() + "): " + sub);
			Global.log().info("freq: " + frequencies.get(i));
			
			double max = Double.NEGATIVE_INFINITY;

			Global.log().info("null model: ER");
			{
				double sizeER = MotifSearchModel.sizeER(data, sub, occs, resets);
				double factorER = baselineER - sizeER;
				factorsER.add(factorER);
				
				max = Math.max(max, factorER);
				 
				Global.log().info("ER baseline: " + baselineER);
				Global.log().info("ER motif code: " + sizeER);
				Global.log().info("ER factor: " + factorER);
			}

			Global.log().info("null model: EL");
			{
				double sizeEL = MotifSearchModel.sizeEL(data, sub, occs,  resets);
					
				double factorEL = baselineEL - sizeEL;
				factorsEL.add(factorEL);
				
				max = Math.max(max, factorEL);
			 
				Global.log().info("EL baseline: " + baselineEL);
				Global.log().info("EL motif code: " + sizeEL);
				Global.log().info("EL factor: " + factorEL);
			}

			Global.log().info("null model: Beta");
			{

				double sizeBeta = MotifSearchModel.sizeBeta(data, sub, occs, resets, betaIterations, betaAlpha, betaSearchDepth);
				double factorBeta = baselineBeta - sizeBeta;
				factorsBeta.add(factorBeta);
			 
				max = Math.max(max, factorBeta);
				
				Global.log().info("Beta baseline: " + baselineBeta);
				Global.log().info("Beta motif code: " + sizeBeta);
				Global.log().info("Beta factor: " + factorBeta);
			}
			
			maxFactors.add(max);
		}
		
		Comparator<Double> comp = Functions.natural();
		Functions.sort(
				factorsBeta, Collections.reverseOrder(comp),
				(List) frequencies,
				(List) factorsER, 
				(List) factorsEL, 
				(List) factorsBeta, 
				(List) subs);
		
		File numbersFile = new File("numbers.csv");
		
		BufferedWriter numbersWriter = new BufferedWriter(new FileWriter(numbersFile));
		for(int i : series(subs.size()))
			numbersWriter.write(frequencies.get(i) + ", " + factorsER.get(i) + ", " + factorsEL.get(i) + ", " + factorsBeta.get(i) + "\n");		
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
		obj.put("directed", directed);
		obj.put("baseline er", baselineER);
		obj.put("baseline el", baselineEL);
		obj.put("baseline beta", baselineBeta);
		Functions.write(obj.toString(), new File( "metadata.json"));
		
		try
		{
			FileIO.python(new File("."), "scripts/plot.py");
		} catch (Exception e)
		{
			Global.log().warning("Failed to run plot script. The script has been copied to the output directory.  (trace:" + e + ")");
		}
	}
}
