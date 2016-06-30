package nl.peterbloem.motive.exec;

import static nl.peterbloem.kit.Functions.log2Choose;
import static nl.peterbloem.kit.Pair.p;
import static nl.peterbloem.kit.Series.series;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

import org.json.JSONObject;
import org.nodes.Graph;
import org.nodes.Graphs;
import org.nodes.Link;
import org.nodes.Node;
import org.nodes.UGraph;
import org.nodes.UNode;
import org.nodes.UTGraph;
import org.nodes.algorithms.Nauty;
import org.nodes.compression.BinomialCompressor;
import org.nodes.data.Data;
import org.nodes.models.DegreeSequenceModel.Prior;
import org.nodes.models.ERSimpleModel;
import org.nodes.models.EdgeListModel;
import org.nodes.random.RandomGraphs;

import nl.peterbloem.kit.FileIO;
import nl.peterbloem.kit.FrequencyModel;
import nl.peterbloem.kit.Functions;
import nl.peterbloem.kit.Global;
import nl.peterbloem.kit.Pair;
import nl.peterbloem.motive.MotifModel;
import nl.peterbloem.motive.MotifSearchModel;
import nl.peterbloem.motive.UPlainMotifExtractor;

public class Synthetic
{
	
	/**
	 *  Size of the random graph in nodes
	 */
	public int n = 5000;
	/**
	 * Size of the random graph in links
	 */
	public int m = 10000;
	
	/**
	 * Size of the injected motif in nodes
	 */
	public int nPrime = 5;
	/**
	 *  Size of the injected motif in links
	 */
	public int mPrime = 6;
	
	/**
	 *  The number of instances to inject (multiple values)
	 */
	public List<Integer> numsInstances = Arrays.asList(0, 10, 100);
	/**
	 *  Maximum degree of instance nodes.
	 */
	public int maxDegree = 5;
	
	/**
	 *  Number of samples to take to find potential motifs
	 */
	public int motifSamples = 5000;
	
	/**
	 *  Number of times to repeat each epxperiment
	 */
	public int runs = 10;

	// * Fields
	private UGraph<String> sub;
	private Map<Pair<UGraph<String>,Integer>, List<Run>> all = 
			new LinkedHashMap<Pair<UGraph<String>,Integer>, List<Run>>();
	
	// * Encountered subgraphs in the correct order
	private List<UGraph<String>> subs;
	
	// * Sums of all factors seen for each subgraph
	private Map<Integer, FrequencyModel<UGraph<String>>> sumFactors = 
			new LinkedHashMap<Integer, FrequencyModel<UGraph<String>>>();
	
	public void main() throws IOException
	{
		subs = Graphs.allIsoConnected(nPrime, "");
		System.out.println(subs.size()  + " " + subs);
		
		// * Sample the subgraph
		sub = RandomGraphs.random(nPrime, mPrime);
		sub = Graphs.blank(sub, "");
		
		for(int numInstances : numsInstances)
			sumFactors.put(numInstances, new FrequencyModel<UGraph<String>>());
		
		// * perform the experiments
		for(int numInstances : numsInstances)
			for(int run : series(runs))
				run(numInstances, run);
		
		sub = (UGraph<String>)Nauty.canonize(sub);
		
		// * sort the subgraphs
		sort();
		
		// * collate the data and print to CSVs
		out();
		
		try
		{
			FileIO.python(new File("."), "scripts/plot.synthetic.py");
		} catch (Exception e)
		{
			throw new RuntimeException("Failed to run plot script. ", e);
		}
	}

	private void sort()
	{
		List<Double> maxFactors = new ArrayList<Double>(subs.size());
		for(UGraph<String> sub : subs)
		{
			List<Double> sums = new ArrayList<Double>(numsInstances.size());
			for(int numInstances : numsInstances)
				sums.add(sumFactors.get(numInstances).frequency(sub));
			maxFactors.add(absMax(sums));
		}

		Comparator<Double> natural = Functions.natural();
		Functions.sort(maxFactors, Collections.reverseOrder(natural), subs);
	}
	
	/**
	 * Returns the element in the list with the greatest absolute magnitude
	 * @param values
	 * @return
	 */
	private double absMax(List<Double> values)
	{
		double maxMag = Double.NEGATIVE_INFINITY;
		double elem = -1.0;
		
		for(double value : values)
			if(Math.abs(value) > maxMag)
			{
				maxMag = Math.abs(value);
				elem = value;
			}
		
		return elem;
	}

	public void run(int numInstances, int run)
	{	
		// * Sample the compressed graph
		UGraph<String> graph = RandomGraphs.random(n - numInstances * (sub.size() - 1), m - numInstances * sub.numLinks());
		
		Global.log().info("Finished sampling subbed graph");
		
		List<UNode<String>> candidates = new ArrayList<UNode<String>>(graph.size());
		for(UNode<String> node : graph.nodes())
			candidates.add(node);
		
		// * remove all candidates with too high degree
		Iterator<UNode<String>> it = candidates.iterator();
		while(it.hasNext())
			if(it.next().degree() > maxDegree)
				it.remove();
		
		if(candidates.size() < numInstances)
			throw new RuntimeException("There are only "+candidates.size()+" potential instance nodes, but "+numInstances+" are required to become instances");
		
		// * Sample a random subset of instances
		List<UNode<String>> instances = Functions.subset(candidates, numInstances);
		
		// * Sample a random multinomial
		List<Double> probs = Functions.randomMultinomial(nPrime);
		
		// * Wire in the instances
		for(Node<String> instance : instances)
		{
			List<Node<String>> newNodes = new ArrayList<Node<String>>();
			// * Add motif nodes
			for(Node<String> subNode : sub.nodes())
				newNodes.add(graph.add(subNode.label()));
			
			// * Wire them up internally
			for(Link<String> link : sub.links())
			{
				int i = link.first().index();
				int j = link.second().index();
				
				newNodes.get(i).connect(newNodes.get(j));
			}
			
			// * Wire them up externally
			for(Link<String> link : instance.links())
			{
				Node<String> other = link.other(instance);
				// - choose a random node inside the motif
				int inside = Functions.choose(probs, 1.0);
				other.connect(newNodes.get(inside));
			}
			
			// * remove the instance node
			instance.remove();
		}
		
		// - Not really necessary...
		graph = Graphs.blank(graph, "");
		
		Global.log().info("Graph reconstructed. nodes: " + graph.size() + ", links: " + graph.numLinks());
		
		// * Perform the motif search
		UPlainMotifExtractor<String> ex = new UPlainMotifExtractor<String>(graph, motifSamples, nPrime);
		
		double baseline = new ERSimpleModel(false).codelength(graph);
		// double baseline = new EdgeListModel(Prior.ML).codelength(graph);
		Global.log().info("baseline " + baseline);
		
//		int nn = graph.size(), nl = graph.numLinks();
//		Global.log().info("choose: " + log2Choose(nl, (nn*nn-nn)/2 ));
		
		for(UGraph<String> s : subs)
		{
			List<List<Integer>> occurrences =  ex.occurrences(s);
			if(occurrences == null)
				occurrences = Collections.emptyList();
			
			Global.log().info("Analysing sub: " + s);
			double motifSize = MotifSearchModel.sizeER(graph, s, occurrences, true);
			Global.log().info("motif size: " + motifSize);
			double factor = baseline - motifSize;
			
			sumFactors.get(numInstances).add(s, factor);
			
			new Run(s, numInstances, run, ex.frequency(s), factor);
		}
	}

	private class Run
	{
		UGraph<String> sub;
		int instances;
		int run;
		
		double frequency;
		double factor;
		
		public Run(UGraph<String> sub, int instances, int run, double frequency, double factor)
		{
			super();
			this.sub = sub;
			this.instances = instances;
			this.run = run;
			
			this.frequency = frequency;
			this.factor = factor;
			
			Pair<UGraph<String>, Integer> pair = Pair.p(sub, instances);
			if(! all.containsKey(pair))
			{
				List<Run> list = new ArrayList<Run>(runs);
				for(int i : series(runs))
					list.add(null);
				all.put(pair, list);
			}
			
			all.get(pair).set(run, this);
		}
		
		public int instances()
		{
			return instances;
		}
		
		public UGraph<String> sub()
		{
			return sub;
		}
		public int run()
		{
			return run;
		}
		
		public double frequency()
		{
			return frequency;
		}
		
		public double factor()
		{
			return factor;
		}

		@Override
		public String toString()
		{
			return  frequency + "_" + factor;
		}
	}

	private void out() throws IOException
	{
		File frequenciesFile = new File("frequencies.csv");
		File factorsFile = new File("factors.csv");
		File meansFile = new File("means.csv");
		
		BufferedWriter frequencies = new BufferedWriter(new FileWriter(frequenciesFile));
		BufferedWriter factors = new BufferedWriter(new FileWriter(factorsFile));
		BufferedWriter means = new BufferedWriter(new FileWriter(meansFile));
		
		// System.out.println(subs);
		
		int i = 0, subIndex = -1;
		for(UGraph<String> sub : subs)
		{
			System.out.println(sub);
			
			if(sub.equals(this.sub))
				subIndex = i;

			// * Write the subgraph
			File graphFile = new File(String.format("motif.%03d.edgelist", i));
			Data.writeEdgeList(sub, graphFile);
			
			// * Write the frequencies, 
			FrequencyModel<Integer> freqSums = new FrequencyModel<Integer>();
			FrequencyModel<Integer> factSums = new FrequencyModel<Integer>();
			
			int c = 0;
			for(int numInstances : numsInstances)
				for(int runIndex : series(runs))
				{
					Run run = all.get(p(sub, numInstances)).get(runIndex);
					frequencies.write((c==0 ? "" : ", ") + run.frequency());
					factors.write((c++==0 ? "" : ", ") + run.factor());
					
					freqSums.add(numInstances, run.frequency());
					factSums.add(numInstances, run.factor());
				}
			
			frequencies.write("\n");
			factors.write("\n");
			
			c = 0;
			for(int numInstances : numsInstances)
				means.write((c++==0 ? "" : ", ") + (freqSums.frequency(numInstances) / (double)runs));
			
			for(int numInstances : numsInstances)
				means.write( ", " + (factSums.frequency(numInstances) / (double)runs));
			means.write("\n");
			
			i++;
		}
		
		JSONObject obj = new JSONObject();
		obj.put("subindex", subIndex);
		obj.put("nums instances", numsInstances);
		obj.put("motif size", nPrime);

		Functions.write(obj.toString(), new File("metadata.json"));
		
		frequencies.close();
		means.close();
		factors.close();
	}
}
