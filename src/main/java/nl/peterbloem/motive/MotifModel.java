package nl.peterbloem.motive;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static nl.peterbloem.kit.Functions.log2;
import static nl.peterbloem.kit.Functions.log2Choose;
import static nl.peterbloem.kit.Functions.log2Factorial;
import static nl.peterbloem.kit.Functions.logFactorial;
import static nl.peterbloem.kit.Functions.max;
import static nl.peterbloem.kit.Functions.prefix;
import static nl.peterbloem.kit.Pair.p;
import static nl.peterbloem.kit.Series.series;
import static org.nodes.motifs.MotifCompressor.MOTIF_SYMBOL;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;

import javax.jws.WebParam.Mode;

import org.nodes.DGraph;
import org.nodes.DLink;
import org.nodes.DNode;
import org.nodes.Graph;
import org.nodes.Graphs;
import org.nodes.Link;
import org.nodes.MapDTGraph;
import org.nodes.MapUTGraph;
import org.nodes.Node;
import org.nodes.UGraph;
import org.nodes.ULink;
import org.nodes.UNode;
import org.nodes.models.DegreeSequenceModel.Prior;
import org.nodes.motifs.MotifCompressor;
import org.nodes.models.DSequenceEstimator;
import org.nodes.models.DegreeSequenceModel;
import org.nodes.models.ERSimpleModel;
import org.nodes.models.EdgeListModel;
import org.nodes.models.Model;
import org.nodes.models.RestrictedToSimple;
import org.nodes.models.StructureModel;
import org.nodes.models.USequenceEstimator;
import static org.nodes.models.DSequenceEstimator.D;
import org.nodes.util.bootstrap.LogNormalCI;

import nl.peterbloem.kit.FrequencyModel;
import nl.peterbloem.kit.Functions;
import nl.peterbloem.kit.Global;
import nl.peterbloem.kit.OnlineModel;
import nl.peterbloem.kit.Pair;
import nl.peterbloem.kit.Series;

/**
 * 
 * Note: we can probably speed up the functions sizeSubbedER and sizeSubbedEL by
 * implementing it with a loop over motif instances, rather than a loop over the 
 * whole graph. This would make the complexity independent of the size of the 
 * graph, allowing graphs of arbitrary size to be tackled.
 * 
 * @author Peter
 *
 * @param <G>
 */
public class MotifModel
{
	private static ExecutorService executor = null;
	
	/**
	 * Sets the threadpool to use (for the beta model). If not set, the beta 
	 * model will manage its own threads.
	 *  
	 * @param executor
	 */
	public static void setExecutor(ExecutorService executor)
	{
		MotifModel.executor = executor;
	}

	
	public static <L> double size(Graph<L> graph, Graph<L> sub,
			List<List<Integer>> occurrences, StructureModel<Graph<?>> nullModel, boolean resetWiring)
	{
		FrequencyModel<String> bits = new FrequencyModel<String>();

		boolean directed = (graph instanceof DGraph<?>); 

		List<List<Integer>> wiring = new ArrayList<List<Integer>>();
		Set<Integer> motifNodes = new HashSet<Integer>();
		Graph<L> subbed;
		if(directed)
			subbed = subbedGraph((DGraph<L>) graph, occurrences, wiring, motifNodes);
		else
			subbed = subbedGraph((UGraph<L>) graph, occurrences, wiring, motifNodes);
		
		if(nullModel instanceof RestrictedToSimple)
		{
			FrequencyModel<Pair<Integer, Integer>> removals = new FrequencyModel<Pair<Integer,Integer>>();

			if(directed)
				subbed = Graphs.toSimpleDGraph((DGraph<L>)subbed, removals);
			else
				subbed = Graphs.toSimpleUGraph((UGraph<L>)subbed, removals);
			
			List<Integer> additions = new ArrayList<Integer>();
			
			for(Link<L> link : subbed.links())
				if(motifNodes.contains(link.first().index()) || motifNodes.contains((link.second().index())))
				{
					int i = link.first().index(), j = link.second().index();
					
					Pair<Integer, Integer> pair = 
							directed ? Pair.p(i, j) : Pair.p(min(i,  j), max(i, j));
				
					additions.add((int)removals.frequency(pair));
				}
						
			bits.add("multiple-edges", Functions.prefix(additions.isEmpty() ? 0 : (long)Functions.max(additions)));
			bits.add("multiple-edges", OnlineModel.storeIntegers(additions)); 
		}
		
		// * Store the labels
		bits.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), subbed.size())); 
								
		bits.add("sub", nullModel.codelength(sub));

		bits.add("subbed", nullModel.codelength(subbed));
		
		// * Store the rewiring information
		bits.add("wiring", wiringBits(sub, wiring, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data 
		bits.add("insertions", log2Factorial(graph.size()) - log2Factorial(subbed.size()));
		
		return bits.total();
	}
	
	public static <L> double sizeBeta(Graph<L> graph, Graph<L> sub,
			List<List<Integer>> occurrences, boolean resetWiring, int iterations, double alpha)
	{		
		if(graph instanceof DGraph<?>)
			return sizeBeta((DGraph<L>)graph, (DGraph<L>) sub, occurrences, resetWiring, iterations, alpha);
		else 
			return sizeBeta((UGraph<L>)graph, (UGraph<L>) sub, occurrences, resetWiring, iterations, alpha);
	}
	
	public static <L> double sizeBeta(DGraph<L> graph, DGraph<L> sub,
			List<List<Integer>> occurrences, boolean resetWiring, int iterations, double alpha)
	{		
		int numThreads = Global.numThreads();
		
		List<List<Integer>> wiring = new ArrayList<List<Integer>>();
		Set<Integer> motifNodes = new HashSet<Integer>();

		// * The "rest" of the code (ie. everything but the parts estimated by 
		//   importance sampling)
		FrequencyModel<String> rest = new FrequencyModel<String>();
		
		// * Compute the degree sequence of the subbed graph
		// ... and any node pairs with multiple links
		List<D> degrees = subbedDegrees(graph, occurrences, rest);
				
		// * The estimated cost of storing the structure of the motif and the 
		//   structure of the subbed graph. 
		List<Double> samples = new ArrayList<Double>(iterations);
		DSequenceEstimator<String> motifModel = new DSequenceEstimator<String>(sub);
		DSequenceEstimator<String> subbedModel = new DSequenceEstimator<String>(degrees);
		motifModel.nonuniform(iterations, numThreads, executor);
		subbedModel.nonuniform(iterations, numThreads, executor);
		
		for(int i : series(iterations))
			samples.add(motifModel.logSamples().get(i) + subbedModel.logSamples().get(i));
		 
		LogNormalCI ci = new LogNormalCI(samples);
				
		// * parameters
		rest.add("sub", DegreeSequenceModel.prior((DGraph<?>)sub, Prior.COMPLETE));
		rest.add("subbed", DegreeSequenceModel.prior(degrees, Prior.COMPLETE));
		
		// * Store the labels
		rest.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), degrees.size())); 
				
		// * Store the rewiring information
		rest.add("wiring", wiringBitsDirect(graph, sub, occurrences, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data 
		rest.add("insertions", log2Factorial(graph.size()) - log2Factorial(degrees.size()));
		
//		System.out.println("ci : " + ci.upperBound(alpha));
//		rest.print(System.out);
		
		return ci.upperBound(alpha) + rest.total();
	}	
	
	public static List<D> subbedDegrees(
			DGraph<?> graph, List<List<Integer>> occurrences, 
			FrequencyModel<String> rest)
	{
		// * records which node is in which occurrence (if any)
		Map<Integer, Integer> nodeInOccurrence = new HashMap<Integer, Integer>();
		
		for(int occurrenceIndex : series(occurrences.size()))
			for(int nodeIndex : occurrences.get(occurrenceIndex))
				nodeInOccurrence.put(nodeIndex, occurrenceIndex);
		
		FrequencyModel<Pair<Integer, Integer>> nodeToInstance = 
				new FrequencyModel<Pair<Integer,Integer>>();
		FrequencyModel<Pair<Integer, Integer>> instanceToNode = 
				new FrequencyModel<Pair<Integer,Integer>>();
		FrequencyModel<Pair<Integer, Integer>> instanceToInstance = 
				new FrequencyModel<Pair<Integer,Integer>>();
		
		// * Record the in and out degrees
		FrequencyModel<Integer> in = new FrequencyModel<Integer>();
		FrequencyModel<Integer> out = new FrequencyModel<Integer>();
		
		for(DLink<?> link : graph.links())
		{
			int fromInstance = nodeInOccurrence.get(link.from().index()) == null ? -1 : nodeInOccurrence.get(link.from().index()); 
			int toInstance =   nodeInOccurrence.get(link.to().index()) == null ? -1 : nodeInOccurrence.get(link.to().index()); 
		
			if(fromInstance == -1 && toInstance == -1)
			{
				out.add(link.from().index());
				in.add(link.to().index());
				continue;
			}
			
			if(fromInstance == -1)
			{
				Pair<Integer, Integer> n2i = Pair.p(link.from().index(), toInstance);
				if(nodeToInstance.frequency(n2i) == 0.0)
				{
					out.add(link.from().index());
					in.add(-(toInstance + 1));
				}
				nodeToInstance.add(n2i);
				continue;
			}
			
			if(toInstance == -1)
			{
				Pair<Integer, Integer> i2n = Pair.p(fromInstance, link.to().index());
				if(instanceToNode.frequency(i2n) == 0.0)
				{
					in.add(link.to().index());
					out.add(-(fromInstance + 1));
				}
				instanceToNode.add(i2n);
				continue;
			}
			
			if(fromInstance != toInstance)
			{
				Pair<Integer, Integer> i2i = Pair.p(fromInstance, toInstance);
				if(instanceToInstance.frequency(i2i) == 0.0)
				{
					out.add(-(fromInstance + 1));
					in.add(-(toInstance + 1));
				}
				instanceToInstance.add(i2i);
			}
		}
		
		Set<Integer> nodes = new LinkedHashSet<Integer>();
		nodes.addAll(in.tokens());
		nodes.addAll(out.tokens());
		List<D> degrees = new ArrayList<D>(nodes.size());
		
		for(int node : nodes)
			degrees.add(new DSequenceEstimator.D((int)in.frequency(node), (int)out.frequency(node)));
				
		List<Integer> additions = new ArrayList<Integer>(graph.size());
		for(Pair<Integer, Integer> token : nodeToInstance.tokens())
			additions.add((int)nodeToInstance.frequency(token) - 1);
		for(Pair<Integer, Integer> token : instanceToNode.tokens())
			additions.add((int)instanceToNode.frequency(token) - 1);
		for(Pair<Integer, Integer> token : instanceToInstance.tokens())
			additions.add((int)instanceToInstance.frequency(token) - 1);
		
		rest.add("multi-edges", Functions.prefix(additions.isEmpty() ? 0 : (long)Functions.max(additions)));
		rest.add("multi-edges", OnlineModel.storeIntegers(additions)); 
		
		// * check for any disconnected nodes and add 0s
		if(! occurrences.isEmpty())
		{	
			int expSize = graph.size() - occurrences.size() * (occurrences.get(0).size() - 1);
			while(degrees.size() < expSize)
				degrees.add(new D(0, 0));
		}
		
		return degrees;
	}
	
	public static <L> double sizeBeta(UGraph<L> graph, UGraph<L> sub,
			List<List<Integer>> occurrences, boolean resetWiring, int iterations, double alpha)
	{		
		int numThreads = Runtime.getRuntime().availableProcessors();
		
		List<List<Integer>> wiring = new ArrayList<List<Integer>>();
		Set<Integer> motifNodes = new HashSet<Integer>();

		// * The "rest" of the code (ie. everything but the parts estimated by 
		//   importance sampling)
		FrequencyModel<String> rest = new FrequencyModel<String>();
		
		// * Compute the degree sequence of the subbed graph
		// ... and any node pairs with multiple links
		List<Integer> degrees = subbedDegrees(graph, occurrences, rest);
				
		// * The estimated cost of storing the structure of the motif and the 
		//   structure of the subbed graph. 
		List<Double> samples = new ArrayList<Double>(iterations);
		USequenceEstimator<String> motifModel = new USequenceEstimator<String>(sub);
		USequenceEstimator<String> subbedModel = new USequenceEstimator<String>(degrees);
		motifModel.nonuniform(iterations, numThreads, executor);
		subbedModel.nonuniform(iterations, numThreads, executor);
		
		for(int i : series(iterations))
			samples.add(motifModel.logSamples().get(i) + subbedModel.logSamples().get(i));
		 
		LogNormalCI ci = new LogNormalCI(samples);
				
		// * parameters
		rest.add("sub", DegreeSequenceModel.prior(sub, Prior.COMPLETE));
		rest.add("subbed", DegreeSequenceModel.priorDegrees(degrees, Prior.COMPLETE));
		
		// * Store the labels
		rest.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), degrees.size())); 
				
		// * Store the rewiring information
		rest.add("wiring", wiringBitsDirect(graph, sub, occurrences, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data 
		rest.add("insertions", log2Factorial(graph.size()) - log2Factorial(degrees.size()));
		
//		System.out.println("ci : " + ci.upperBound(alpha));
//		rest.print(System.out);
		
		return ci.upperBound(alpha) + rest.total();
	}		
	
	public static List<Integer> subbedDegrees(
			UGraph<?> graph, List<List<Integer>> occurrences, 
			FrequencyModel<String> rest)
	{
		// * records which node is in which occurrence (if any)
		Map<Integer, Integer> nodeInOccurrence = new HashMap<Integer, Integer>();
		
		for(int occurrenceIndex : series(occurrences.size()))
			for(int nodeIndex : occurrences.get(occurrenceIndex))
				nodeInOccurrence.put(nodeIndex, occurrenceIndex);
		
		FrequencyModel<Pair<Integer, Integer>> nodeToInstance = 
				new FrequencyModel<Pair<Integer,Integer>>();
		FrequencyModel<Pair<Integer, Integer>> instanceToInstance = 
				new FrequencyModel<Pair<Integer,Integer>>();
		
		// * Record the in and out degrees
		FrequencyModel<Integer> degrees = new FrequencyModel<Integer>();

		
		for(ULink<?> link : graph.links())
		{
			int firstInstance = nodeInOccurrence.get(link.first().index()) == null ? -1 : nodeInOccurrence.get(link.first().index()); 
			int secondInstance =   nodeInOccurrence.get(link.second().index()) == null ? -1 : nodeInOccurrence.get(link.second().index()); 
		
			if(firstInstance == -1 && secondInstance == -1)
			{
				degrees.add(link.first().index());
				degrees.add(link.second().index());
				continue;
			}
			
			if(firstInstance == -1)
			{
				Pair<Integer, Integer> n2i = Pair.p(link.first().index(), secondInstance);
				if(nodeToInstance.frequency(n2i) == 0.0)
				{
					degrees.add(link.first().index());
					degrees.add(-(secondInstance + 1));
				}
				nodeToInstance.add(n2i);
				continue;
			}
			
			if(secondInstance == -1)
			{
				Pair<Integer, Integer> n2i = Pair.p(link.second().index(), firstInstance);
				if(nodeToInstance.frequency(n2i) == 0.0)
				{
					degrees.add(link.second().index());
					degrees.add(-(firstInstance + 1));
				}
				nodeToInstance.add(n2i);
				continue;
			}
			
			if(firstInstance != secondInstance)
			{
				Pair<Integer, Integer> i2i = Pair.p(min(firstInstance, secondInstance), max(firstInstance, secondInstance));
				if(instanceToInstance.frequency(i2i) == 0.0)
				{
					degrees.add(-(firstInstance + 1));
					degrees.add(-(secondInstance + 1));
				}
				instanceToInstance.add(i2i);
			}
		}
		
		List<Integer> result = new ArrayList<Integer>((int)degrees.distinct());
		
		for(int node : degrees.tokens())
			result.add((int)degrees.frequency(node));
				
		List<Integer> additions = new ArrayList<Integer>(graph.size());
		for(Pair<Integer, Integer> token : nodeToInstance.tokens())
			additions.add((int)nodeToInstance.frequency(token) - 1);
		for(Pair<Integer, Integer> token : instanceToInstance.tokens())
			additions.add((int)instanceToInstance.frequency(token) - 1);
		
		rest.add("multi-edges", Functions.prefix(
				additions.isEmpty() ? 
				0 : (long)Functions.max(additions)));
		rest.add("multi-edges", OnlineModel.storeIntegers(additions)); 
		
		// * check for any disconnected nodes and add 0s
		if(! occurrences.isEmpty())
		{
			int expSize = graph.size() - occurrences.size() * (occurrences.get(0).size() - 1);
			while(result.size() < expSize)
				result.add(0);
		}
		
		return result;
	}	
	
	public static double wiringBits(Graph<?> sub, List<List<Integer>> wiring,
			boolean reset)
	{
		OnlineModel<Integer> om = new OnlineModel<Integer>(Series.series(sub.size()));

		double bits = 0.0;
		for (List<Integer> motifWires : wiring)
		{
			if (reset)
				om = new OnlineModel<Integer>(Series.series(sub.size()));

			for (int wire : motifWires)
				bits += - Functions.log2(om.observe(wire));
		}
		
		return bits;
	}

	public static double sizeER(Graph<?> graph, Graph<?> sub, List<List<Integer>> occurrences, boolean resetWiring)
	{
		if(graph instanceof DGraph)
			return sizeER((DGraph<?>) graph, (DGraph<?>) sub, occurrences, resetWiring);
		else
			return sizeER((UGraph<?>) graph, (UGraph<?>) sub, occurrences, resetWiring); 
	}
	
	private static ERSimpleModel erModel = new ERSimpleModel(true);
	
	public static double sizeER(DGraph<?> graph, DGraph<?> sub,
				List<List<Integer>> occurrences, boolean resetWiring)
	{		
		FrequencyModel<String> bits = new FrequencyModel<String>();
		
		bits.add("sub", erModel.codelength(sub));
		sizeSubbedER(graph, sub, occurrences, bits);
		
		// * Store the rewiring information
		bits.add("wiring", wiringBitsDirect(graph, sub, occurrences, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data
		int subbedSize = graph.size() - (sub.size() - 1) * occurrences.size(); 
		
		bits.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), subbedSize)); 

		bits.add("insertions", log2Factorial(graph.size()) - log2Factorial(subbedSize));
		
		return bits.total();
	}
	
	public static double sizeER(UGraph<?> graph, UGraph<?> sub,
				List<List<Integer>> occurrences, boolean resetWiring)
	{		
		FrequencyModel<String> bits = new FrequencyModel<String>();
		
		bits.add("sub", erModel.codelength(sub));
		sizeSubbedER(graph, sub, occurrences, bits);
		
		// * Store the rewiring information
		bits.add("wiring", wiringBitsDirect(graph, sub, occurrences, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data
		int subbedSize = graph.size() - (sub.size() - 1) * occurrences.size(); 
		bits.add("insertions", log2Factorial(graph.size()) - log2Factorial(subbedSize));
		bits.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), subbedSize)); 
		
		// bits.print(System.out);
		
		return bits.total();
	}
	
	/**
	 * A version of the ER model that loops only over the instances. It requires 
	 * the degrees of the graph to be given.
	 * 
	 * @param graph
	 * @param degrees
	 * @param sub
	 * @param occurrences
	 * @param resetWiring
	 * @return
	 */
	public static double sizeERInst(UGraph<?> graph, UGraph<?> sub,
			List<List<Integer>> occurrences, boolean resetWiring)
	{		
		FrequencyModel<String> bits = new FrequencyModel<String>();
		
		bits.add("sub", erModel.codelength(sub));

		FrequencyModel<Pair<Integer, Integer>> multiEdges = new FrequencyModel<Pair<Integer, Integer>>();
		List<List<Integer>> rewiring = new LinkedList<List<Integer>>();
		
		Pair<Long, Long> pair = subbedERInstances(graph, sub, occurrences, multiEdges, rewiring);
		
		// * store the template graph (as a simple graph) 
		bits.add("subbed", ERSimpleModel.undirected(pair.first(), pair.second(), true));
		
		// * store the multi-edges
		// - We are storing, for each link, the number of additional edges required
		//   (so everything's - 1)
		bits.add("multi-edges", multiEdges(multiEdges));
		
		// * Store the rewiring information
		bits.add("wiring", wiringBits(sub, rewiring, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data
		int subbedSize = graph.size() - (sub.size() - 1) * occurrences.size();
				
		bits.add("insertions", log2Factorial(graph.size()) - log2Factorial(subbedSize));
		bits.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), subbedSize)); 
				
		// bits.print(System.out);
		
		return bits.total();
	}
	
	/**
	 * A version of the EL model that loops only over the instances. It requires 
	 * the degrees of the graph to be given.
	 * 
	 * @param graph
	 * @param degrees
	 * @param sub
	 * @param occurrences
	 * @param resetWiring
	 * @return
	 */
	public static double sizeERInst(DGraph<?> graph, DGraph<?> sub,
			List<List<Integer>> occurrences, boolean resetWiring)
	{		
		FrequencyModel<String> bits = new FrequencyModel<String>();
		
		bits.add("sub", erModel.codelength(sub));

		FrequencyModel<Pair<Integer, Integer>> multiEdges = new FrequencyModel<Pair<Integer, Integer>>();
		List<List<Integer>> rewiring = new LinkedList<List<Integer>>();
		
		Pair<Long, Long> pair = subbedERInstances(graph, sub, occurrences, multiEdges, rewiring);
		
		// * store the template graph (as a simple graph) -
		bits.add("subbed", ERSimpleModel.directed(pair.first(), pair.second(), true));
		
		// * store the multi-edges
		// - We are storing, for each link, the number of additional edges required
		//   (so everything's - 1)
		bits.add("multi-edges", multiEdges(multiEdges));
		
		// * Store the rewiring information
		bits.add("wiring", wiringBits(sub, rewiring, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data
		long subbedSize = ((long)graph.size()) - (sub.size() - 1) * (long)occurrences.size();
				
		bits.add("insertions", log2Factorial(graph.size()) - log2Factorial(subbedSize));
		bits.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), subbedSize)); 
				
		// bits.print(System.out);
		
		return bits.total();
	}	
	
	private static void sizeSubbedER(DGraph<?> graph, DGraph<?> sub,
			List<List<Integer>> occurrences, FrequencyModel<String> bits)
	{
		int subbedSize = graph.size() - (sub.size() - 1) * occurrences.size();
		int subbedLinks = 0;
		// - NB: subbedLinks is not simply 
		//      graph.numLinks() - (sub.numLinks() * occurrences.size())
		//   because the subbed graph is a simple graph: any multiple edges 
		//   created by removal of occurrences are removed and stored elsewhere. 
		
		// * records which node is in which occurrence (if any)
		Map<Integer, Integer> nodeInOccurrence = new HashMap<Integer, Integer>();
		
		for(int occurrenceIndex : series(occurrences.size()))
			for(int nodeIndex : occurrences.get(occurrenceIndex))
				nodeInOccurrence.put(nodeIndex, occurrenceIndex);
		
		FrequencyModel<Pair<Integer, Integer>> nodeToInstance = 
				new FrequencyModel<Pair<Integer,Integer>>();
		FrequencyModel<Pair<Integer, Integer>> instanceToNode = 
				new FrequencyModel<Pair<Integer,Integer>>();
		FrequencyModel<Pair<Integer, Integer>> instanceToInstance = 
				new FrequencyModel<Pair<Integer,Integer>>();
		
		for(DLink<?> link : graph.links())
		{
			int fromInstance = nodeInOccurrence.get(link.from().index()) == null ? -1 : nodeInOccurrence.get(link.from().index()); 
			int toInstance =   nodeInOccurrence.get(link.to().index()) == null ? -1 : nodeInOccurrence.get(link.to().index()); 
		
			if(fromInstance == -1 && toInstance == -1)
			{
				subbedLinks++;
				continue;
			}
			
			if(fromInstance == -1)
			{
				Pair<Integer, Integer> n2i = Pair.p(link.from().index(), toInstance);
				if(nodeToInstance.frequency(n2i) == 0.0)
					subbedLinks++;
				nodeToInstance.add(n2i);
				continue;
			}
			
			if(toInstance == -1)
			{
				Pair<Integer, Integer> i2n = Pair.p(fromInstance, link.to().index());
				if(instanceToNode.frequency(i2n) == 0.0)
					subbedLinks++;
				instanceToNode.add(i2n);
				continue;
			}
			
			if(fromInstance != toInstance)
			{
				Pair<Integer, Integer> i2i = Pair.p(fromInstance, toInstance);
				if(instanceToInstance.frequency(i2i) == 0.0)
					subbedLinks++;
				instanceToInstance.add(i2i);
			}
		}
		
		// * size of the subbed graph under the binomial compressor
		bits.add("subbed", ERSimpleModel.directed(subbedSize, subbedLinks, true));
		
		List<Integer> additions = new ArrayList<Integer>(graph.size());
		for(Pair<Integer, Integer> token : nodeToInstance.tokens())
			additions.add((int)nodeToInstance.frequency(token) - 1);
		for(Pair<Integer, Integer> token : instanceToNode.tokens())
			additions.add((int)instanceToNode.frequency(token) - 1);
		for(Pair<Integer, Integer> token : instanceToInstance.tokens())
			additions.add((int)instanceToInstance.frequency(token) - 1);
		
		bits.add("multiple-edges", Functions.prefix(additions.isEmpty() ? 0 : (long)max(additions)));
		bits.add("multiple-edges", OnlineModel.storeIntegers(additions)); 
	}
	
	private static void sizeSubbedER(UGraph<?> graph, UGraph<?> sub,
			List<List<Integer>> occurrences, FrequencyModel<String> bits)
	{
		int subbedSize = graph.size() - (sub.size() - 1) * occurrences.size();
		int subbedLinks = 0;
		
		// * records which node is in which occurrence (if any)
		Map<Integer, Integer> nodeInOccurrence = new HashMap<Integer, Integer>();
		
		for(int occurrenceIndex : series(occurrences.size()))
			for(int nodeIndex : occurrences.get(occurrenceIndex))
				nodeInOccurrence.put(nodeIndex, occurrenceIndex);
		
		FrequencyModel<Pair<Integer, Integer>> nodeToInstance = 
				new FrequencyModel<Pair<Integer,Integer>>();
		FrequencyModel<Pair<Integer, Integer>> instanceToInstance = 
				new FrequencyModel<Pair<Integer,Integer>>();
		
		for(Link<?> link : graph.links())
		{
			int firstInstance = nodeInOccurrence.get(link.first().index()) == null ? -1 : nodeInOccurrence.get(link.first().index()); 
			int secondInstance =   nodeInOccurrence.get(link.second().index()) == null ? -1 : nodeInOccurrence.get(link.second().index()); 
		
			if(firstInstance == -1 && secondInstance == -1)
			{
				subbedLinks++;
				continue;
			}
			
			if(firstInstance == -1) // second is in an instance
			{
				Pair<Integer, Integer> n2i = Pair.p(link.first().index(), secondInstance);
				if(nodeToInstance.frequency(n2i) == 0.0)
					subbedLinks++;
				nodeToInstance.add(n2i);
				continue;
			}
			
			if(secondInstance == -1) // first is in an instance
			{
				Pair<Integer, Integer> n2i = Pair.p(link.second().index(), firstInstance);
				if(nodeToInstance.frequency(n2i) == 0.0)
					subbedLinks++;
				nodeToInstance.add(n2i);
				continue;
			}
			
			if(firstInstance != secondInstance) // both are in an instance
			{
				Pair<Integer, Integer> i2i = Pair.p(Math.min(firstInstance, secondInstance), max(firstInstance, secondInstance));
				if(instanceToInstance.frequency(i2i) == 0.0)
					subbedLinks++;
				instanceToInstance.add(i2i);
			}
		}
		
		// * size of the subbed graph under the binomial compressor
		bits.add("subbed", ERSimpleModel.undirected(subbedSize, subbedLinks, true));
				
		List<Integer> additions = new ArrayList<Integer>(graph.size());
		for(Pair<Integer, Integer> token : nodeToInstance.tokens())
			additions.add((int)nodeToInstance.frequency(token) - 1);
		for(Pair<Integer, Integer> token : instanceToInstance.tokens())
			additions.add((int)instanceToInstance.frequency(token) - 1);
				
		bits.add("multiple-edges", Functions.prefix(additions.isEmpty() ? 0 : (long)max(additions)));
		bits.add("multiple-edges", OnlineModel.storeIntegers(additions)); 
	}	

	public static double sizeEL(Graph<?> graph, Graph<?> sub, List<List<Integer>> occurrences, boolean resetWiring)
	{
		if(graph instanceof DGraph)
			return sizeEL((DGraph<?>) graph, (DGraph<?>) sub, occurrences, resetWiring);
		else
			return sizeEL((UGraph<?>) graph, (UGraph<?>) sub, occurrences, resetWiring); 
	}
	
	private static EdgeListModel elModel = new EdgeListModel(Prior.COMPLETE);
	public static double sizeEL(DGraph<?> graph, DGraph<?> sub,
			List<List<Integer>> occurrences, boolean resetWiring)
	{		
		FrequencyModel<String> bits = new FrequencyModel<String>();
		
		bits.add("sub", elModel.codelength(sub));

		List<D> degrees = subbedDegrees(graph, occurrences, bits);
		bits.add("subbed", EdgeListModel.directed(degrees, Prior.COMPLETE));
		
		// * Store the rewiring information
		bits.add("wiring", wiringBitsDirect(graph, sub, occurrences, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data
		int subbedSize = graph.size() - (sub.size() - 1) * occurrences.size(); 
		bits.add("insertions", log2Factorial(graph.size()) - log2Factorial(subbedSize));
		bits.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), subbedSize)); 
		
		// bits.print(System.out);
		
		return bits.total();
	}
	
	public static double sizeEL(UGraph<?> graph, UGraph<?> sub,
			List<List<Integer>> occurrences, boolean resetWiring)
	{		
		FrequencyModel<String> bits = new FrequencyModel<String>();
		
		bits.add("sub", elModel.codelength(sub));

		List<Integer> degrees = subbedDegrees(graph, occurrences, bits);
		bits.add("subbed", EdgeListModel.undirected(degrees, Prior.COMPLETE));
		
		// * Store the rewiring information
		bits.add("wiring", wiringBitsDirect(graph, sub, occurrences, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data
		int subbedSize = graph.size() - (sub.size() - 1) * occurrences.size(); 
		bits.add("insertions", log2Factorial(graph.size()) - log2Factorial(subbedSize));
		bits.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), subbedSize)); 
				
		// bits.print(System.out);

		return bits.total();
	}
	
	/**
	 * A version of the EL model that loops only over the instances. It requires 
	 * the degrees of the graph to be given.
	 * 
	 * @param graph
	 * @param degrees
	 * @param sub
	 * @param occurrences
	 * @param resetWiring
	 * @return
	 */
	public static double sizeEL(DGraph<?> graph, List<D> degrees, DGraph<?> sub,
			List<List<Integer>> occurrences, boolean resetWiring)
	{		
		FrequencyModel<String> bits = new FrequencyModel<String>();
		
		bits.add("sub", elModel.codelength(sub));

		FrequencyModel<Pair<Integer, Integer>> multiEdges = new FrequencyModel<Pair<Integer, Integer>>();
		List<List<Integer>> rewiring = new LinkedList<List<Integer>>();
		
		List<D> sDegrees = subbedDegrees(graph, degrees, occurrences, multiEdges, rewiring);
		
		// * store the template graph (as a simple graph) 
		bits.add("subbed", EdgeListModel.directed(sDegrees, Prior.COMPLETE));
		
		// * store the multi-edges
		bits.add("multi-edges", multiEdges(multiEdges));
		
		// * Store the rewiring information
		bits.add("wiring", wiringBits(sub, rewiring, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data
		long subbedSize = (long)graph.size() - (sub.size() - 1) * (long)occurrences.size();
		
		assert(sDegrees.size() == subbedSize);
		
		bits.add("insertions", log2Factorial(graph.size()) - log2Factorial(subbedSize));
		bits.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), subbedSize)); 
				
		// bits.print(System.out);
		
		return bits.total();
	}	
	
	public static double multiEdges(FrequencyModel<Pair<Integer, Integer>> multiEdges)
	{
		// - We are storing, for each link, the number of _additional_ edges required
		//   (so everything's - 1)
	
		if(multiEdges.tokens().isEmpty())
			return Functions.prefix(0);
		
		double mBits = 0.0;
		int max = (int)multiEdges.frequency(multiEdges.maxToken());

		mBits += Functions.prefix(max - 1);
		OnlineModel<Integer> model = new OnlineModel<Integer>(Series.series(0, max));

		// - loop over all rewired edges
		for(Pair<Integer, Integer> token : multiEdges.tokens())
			mBits += model.encode((int)multiEdges.frequency(token) - 1);
		
		return mBits;
	}
	
	/**
	 * A version of the EL model that loops only over the instances. It requires 
	 * the degrees of the graph to be given.
	 * 
	 * @param graph
	 * @param degrees
	 * @param sub
	 * @param occurrences
	 * @param resetWiring
	 * @return
	 */
	public static double sizeEL(UGraph<?> graph, List<Integer> degrees, UGraph<?> sub,
			List<List<Integer>> occurrences, boolean resetWiring)
	{		
		FrequencyModel<String> bits = new FrequencyModel<String>();
		
		bits.add("sub", elModel.codelength(sub));

		FrequencyModel<Pair<Integer, Integer>> multiEdges = new FrequencyModel<Pair<Integer, Integer>>();
		List<List<Integer>> rewiring = new LinkedList<List<Integer>>();
		
		List<Integer> sDegrees = subbedDegrees(graph, degrees, occurrences, multiEdges, rewiring);
		
		// * store the template graph (as a simple graph) 
		bits.add("subbed", EdgeListModel.undirected(sDegrees, Prior.COMPLETE));
		
		// * store the multi-edges
		bits.add("multi-edges", multiEdges(multiEdges));
		
		// * Store the rewiring information
		bits.add("wiring", wiringBits(sub, rewiring, resetWiring));
		
		// * Store the insertion order, to preserve the precise ordering of the
		//   nodes in the data
		int subbedSize = graph.size() - (sub.size() - 1) * occurrences.size();
		
		assert(sDegrees.size() == subbedSize);
		
		bits.add("insertions", log2Factorial(graph.size()) - log2Factorial(subbedSize));
		bits.add("labels", Functions.prefix(occurrences.size()) + log2Choose(occurrences.size(), subbedSize)); 
				
		// bits.print(System.out);
		
		return bits.total();
	}

	public static <L> double wiringBitsDirect(Graph<L> graph, Graph<?> sub, List<List<Integer>> occurrences,
			boolean reset)
	{
		OnlineModel<Integer> om = new OnlineModel<Integer>(Series.series(sub.size()));
		
		double wiringBits = 0.0;
		
		for (List<Integer> occurrence : occurrences)
		{			
			if(reset)
				om = new OnlineModel<Integer>(Series.series(sub.size()));
			
			// * The index of the node within the occurrence
			for (int indexInSubgraph : series(occurrence.size()))
			{
				Node<L> node = graph.get(occurrence.get(indexInSubgraph));
				
				for(Link<L> link : node.links())
				{
					Node<L> neighbor = link.other(node);
					
					if(! occurrence.contains(neighbor.index()))
						wiringBits += - log2(om.observe(indexInSubgraph));
				}
			}
		}	
		
		return wiringBits;
	}
	
	/**
	 * Create a copy of the input graph with all (known) occurrences of the 
	 * given subgraph replaced by a single node.
	 * 
	 * @param inputGraph
	 * @param sub
	 * @param occurrences non-overlapping
	 * @param wiring
	 * @return
	 */
	public static <L> DGraph<L> subbedGraph(
			DGraph<L> inputGraph,
			List<List<Integer>> occurrences,
			List<List<Integer>> wiring, Set<Integer> motifNodes)
	{
		// * Create a copy of the input.
		//   We will re-purpose node 0 of each occurrence as the new instance 
		//   node, so for those nodes, we set the label to null
		DGraph<L> copy = new MapDTGraph<L, String>();
		Set<Node<L>> mNodes = new LinkedHashSet<Node<L>>();
		
		Set<Integer> firstNodes = new HashSet<Integer>();
		for(List<Integer> occurrence : occurrences)
			firstNodes.add(occurrence.get(0));
		
		// -- copy the nodes
		for (DNode<L> node : inputGraph.nodes())
			if(firstNodes.contains(node.index()))
				mNodes.add(copy.add(null));
			else
				copy.add(node.label());
		
		// note: slightly leaky abstraction, we are counting on MapDTGraph's 
		// persistent node objects (ie. the Node objects in motif Nodes will
		// stay up to date even if we remove others from the graph later.)
		
		// -- copy the links
		for (DLink<L> link : inputGraph.links())
		{
			int i = link.from().index(), 
			    j = link.to().index();
			
			copy.get(i).connect(copy.get(j));
		}

		// * Translate the occurrences from integers to nodes (in the copy)
		List<List<DNode<L>>> occ = 
				new ArrayList<List<DNode<L>>>(occurrences.size());
		
		for (List<Integer> occurrence : occurrences)
		{
			List<DNode<L>> nodes = 
					new ArrayList<DNode<L>>(occurrence.size());
			for (int index : occurrence)
				nodes.add(copy.get(index));
			
			occ.add(nodes);
		}
		
		for (List<DNode<L>> occurrence : occ)
		{
			// * Use the first node of the motif as the symbol node
			DNode<L> newNode = occurrence.get(0);

			// - This will hold the information how each edge into the motif node should be wired
			//   into the motif subgraph (to be encoded later)
			List<Integer> motifWiring = new ArrayList<Integer>(occ.size());
			wiring.add(motifWiring);
			
			for (int indexInSubgraph : series(occurrence.size()))
			{
				DNode<L> node = occurrence.get(indexInSubgraph);
				
				for(DLink<L> link : node.links())
				{
					// If the link is external
					DNode<L> neighbor = link.other(node);
					
					if(! occurrence.contains(neighbor))
					{
						if(!node.equals(newNode))
						{
							if(link.from().equals(node))
								newNode.connect(neighbor);
							else
								neighbor.connect(newNode);
						}
					
						motifWiring.add(indexInSubgraph);
					}
				}
			}

			for (int i : series(1, occurrence.size()))
				occurrence.get(i).remove();
		}
		
		for(Node<L> node : mNodes)
			motifNodes.add(node.index());

		return copy;
	}
	
	/**
	 * Create a copy of the input graph with all given occurrences of the 
	 * given subgraph replaced by a single node.
	 * 
	 * @param inputGraph
	 * @param sub
	 * @param occurrences
	 * @param wiring
	 * @return
	 */
	public static <L> UGraph<L> subbedGraph(
			UGraph<L> inputGraph,
			List<List<Integer>> occurrences,
			List<List<Integer>> wiring,
			Set<Integer> motifNodes)
	{	
		// * Create a copy of the input.
		//   We will re-purpose node 0 of each occurrence as the new instance 
		//   node, so for those nodes, we set the label to null
		UGraph<L> copy = new MapUTGraph<L, String>();
		Set<Node<L>> mNodes = new LinkedHashSet<Node<L>>();
		
		Set<Integer> firstNodes = new HashSet<Integer>();
		for(List<Integer> occurrence : occurrences)
			firstNodes.add(occurrence.get(0));
		
		// -- copy the nodes
		for (UNode<L> node : inputGraph.nodes())
			if(firstNodes.contains(node.index()))
				mNodes.add(copy.add(null));
			else
				copy.add(node.label());
		
		// -- copy the links
		for (ULink<L> link : inputGraph.links())
		{
			int i = link.first().index(), 
			    j = link.second().index();
			
			copy.get(i).connect(copy.get(j));
		}
		
		// * Translate the occurrences from integers to nodes (in the copy)
		List<List<UNode<L>>> occ = 
				new ArrayList<List<UNode<L>>>(occurrences.size());
		
		for (List<Integer> occurrence : occurrences)
		{
			List<UNode<L>> nodes = 
					new ArrayList<UNode<L>>(occurrence.size());
			for (int index : occurrence)
				nodes.add(copy.get(index));
			
			occ.add(nodes);
		}
		
		for (List<UNode<L>> occurrence : occ)
		{
				// * Wire a new symbol node into the graph to represent the occurrence
				UNode<L> newNode = occurrence.get(0);

				// - This will hold the information how each edge into the motif node should be wired
				//   into the motif subgraph (to be encoded later)
				List<Integer> motifWiring = new ArrayList<Integer>();
				wiring.add(motifWiring);
				
				for (int indexInSubgraph : series(occurrence.size()))
				{
					UNode<L> node = occurrence.get(indexInSubgraph);
					
					for(ULink<L> link : node.links())
					{
						UNode<L> neighbor = link.other(node);
						
						if(! occurrence.contains(neighbor))	// If the link is external
						{
							if(!node.equals(newNode))
								newNode.connect(neighbor);

							motifWiring.add(indexInSubgraph);
						}
					}
				}

				for (int i : series(1, occurrence.size()))
					occurrence.get(i).remove();
		}
		
		for(Node<L> node : mNodes)
			motifNodes.add(node.index());

		return copy;
	}
	
	/**
	 * Computes the degree sequence of the template graph. This method loops only 
	 * over the list of occurrences, making it faster for large graphs with 
	 * few occurrences.
	 * 
	 * @param graph
	 * @param degrees
	 * @param sub
	 * @param occurrences
	 * @param multiEdges An empty frequencymodel receiving how often certain edges in 
	 * the template graph should be repeated (one occurrence in the fm no repeats). 
	 * For performance reasons, the actual indices refer to the old graph, not 
	 * the template graph. 
	 * @param rewiring An empty list, receiving the sequence of rewiring integers.
	 * @return
	 */
	public static List<Integer> subbedDegrees(
			UGraph<?> graph, List<Integer> degrees, 
			List<List<Integer>> occurrences,
			FrequencyModel<Pair<Integer, Integer>> multiEdges,
			List<List<Integer>> rewiring)
	{
		if(occurrences.isEmpty())
			return degrees;

		List<Integer> subbedDegrees = new ArrayList<Integer>(degrees);
		
		// * Which nodes have been mapped to which instance node 
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		// * Links that have been rewired. We store the old link (ie. the left 
		//   and right indices in the original graph). 
		Set<Pair<Integer, Integer>> rewLinks = 
				new LinkedHashSet<Pair<Integer, Integer>>();
		
		for(List<Integer> occurrence : occurrences)
		{	

			List<Integer> rw = new LinkedList<Integer>(); 
			// * Remove the occurrence (first node becomes instance node)
			for(int index : occurrence.subList(1, occurrence.size()))
				subbedDegrees.set(index, -1);
			subbedDegrees.set(occurrence.get(0), 0);

			for(int index : occurrence)
				map.put(index, occurrence.get(0));
			
			// * Remove all links linking into an occurrence
			for(int i : series(occurrence.size()))
			{
				int index = occurrence.get(i);
				for(UNode<?> node : graph.get(index).neighbors())
					if(! occurrence.contains(node.index()))
					{		
						subbedDegrees.set(node.index(), 
							subbedDegrees.get(node.index()) - 1);
						
						rewLinks.add(ordered(index, node.index()));
						
						rw.add(i);						
					}
			}
				
			rewiring.add(rw);
		}
		
		// * set occurrence nodes back to 0
		for(List<Integer> occurrence : occurrences)
			subbedDegrees.set(occurrence.get(0), 0);
				
		// * convert the rewritten links to new indices, and build a 
		//  frequencymodel 
		for(Pair<Integer, Integer> link : rewLinks)
		{
			int f = link.first(), s = link.second();
			int a = map.containsKey(f) ? map.get(f) : f;
			int b = map.containsKey(s) ? map.get(s) : s;
			
			multiEdges.add(ordered(a, b));
		}
		
		// * Add each rewritten link _once_
		for(Pair<Integer, Integer> link : multiEdges.tokens())
		{
			subbedDegrees.set(link.first(),  subbedDegrees.get(link.first())  + 1);
			subbedDegrees.set(link.second(), subbedDegrees.get(link.second()) + 1);
		}
						
		List<Integer> res = new ArrayList<Integer>(
			graph.size() - occurrences.size() * (occurrences.get(0).size() - 1));
		for(int degree : subbedDegrees)
			if(degree >= 0)
				res.add(degree);

		return res;
	}
	
	/**
	 * Computes the degree sequence of the template graph. This method loops only 
	 * over the list of occurrences, making it faster for large graphs with 
	 * few occurrences.
	 * 
	 * @param graph
	 * @param degrees
	 * @param sub
	 * @param occurrences
	 * @param multiEdges An empty frequencymodel receiving how often certain edges in 
	 * the template graph should be repeated (one occurrence in the fm no repeats). 
	 * For performance reasons, the actual indices refer to the old graph, not 
	 * the template graph. 
	 * @param rewiring An empty list, receiving the sequence of rewiring integers.
	 * @return
	 */
	public static List<D> subbedDegrees(
			DGraph<?> graph, List<D> degrees, 
			List<List<Integer>> occurrences,
			FrequencyModel<Pair<Integer, Integer>> multiEdges,
			List<List<Integer>> rewiring)
	{
		if(occurrences.isEmpty())
			return degrees;
	
		List<D> subbedDegrees = new ArrayList<D>(degrees.size());
		for(D degree : degrees)
			subbedDegrees.add(new D(degree.in(), degree.out()));
		
		// * Which nodes have been mapped to which instance node 
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		// * Links that have been rewired. We store the old link (ie. the left 
		//   and right indices in the original graph). 
		Set<Pair<Integer, Integer>> rewLinks = 
				new LinkedHashSet<Pair<Integer, Integer>>();
		
		for(List<Integer> occurrence : occurrences)
		{	
	
			List<Integer> rw = new LinkedList<Integer>(); 
			// * Remove the occurrence (first node becomes instance node)
			for(int index : occurrence.subList(1, occurrence.size()))
				subbedDegrees.set(index, null);
			subbedDegrees.set(occurrence.get(0), new D(0, 0));
	
			for(int index : occurrence)
				map.put(index, occurrence.get(0));
			
			// * Remove all links linking into an occurrence
			for(int i : series(occurrence.size()))
			{
				int index = occurrence.get(i);
				for(DNode<?> node : graph.get(index).out())
					if(! occurrence.contains(node.index()))
					{		
						D old = subbedDegrees.get(node.index());
						if(old != null)
							subbedDegrees.set(node.index(), new D(old.in() - 1, old.out()));
						
						rewLinks.add(Pair.p(index, node.index()));
						
						rw.add(i);						
					}
				
				for(DNode<?> node : graph.get(index).in())
					if(! occurrence.contains(node.index()))
					{		
						D old = subbedDegrees.get(node.index());
						if(old != null)
							subbedDegrees.set(node.index(), new D(old.in(), old.out() - 1));
							
						rewLinks.add(Pair.p(node.index(), index));
						
						rw.add(i);						
					}
			}
				
			rewiring.add(rw);
		}
		
		// * set occurrence nodes back to 0
		for(List<Integer> occurrence : occurrences)
			subbedDegrees.set(occurrence.get(0), new D(0, 0));
				
		// * convert the rewritten links to new indices, and build a 
		//  frequencymodel 
		int size = 0;
		for(Pair<Integer, Integer> link : rewLinks)
		{
			int f = link.first(), s = link.second();
			int a = map.containsKey(f) ? map.get(f) : f;
			int b = map.containsKey(s) ? map.get(s) : s;
			
			multiEdges.add(Pair.p(a, b));
			
			size ++;
			if(size % 10000 == 0)
				System.out.println(size + " rewritten links processed");
		}
	
		System.out.println(".");
		
		// * Add each rewritten link _once_
		for(Pair<Integer, Integer> link : multiEdges.tokens())
		{
			D old;
			old = subbedDegrees.get(link.first());
			subbedDegrees.set(link.first(), new D(old.in(), old.out() + 1));
			
			old = subbedDegrees.get(link.second());
			subbedDegrees.set(link.second(), new D(old.in() + 1, old.out()));
		}
						
		List<D> res = new ArrayList<D>(
			graph.size() - occurrences.size() * (occurrences.get(0).size() - 1));
		
		for(D degree : subbedDegrees)
			if(degree != null)
				res.add(degree);
	
	
		return res;
	}


	/**
	 * Computes the size and number of links in the template graph, by looping 
	 * only over the instances. This method should be fast for large graphs with
	 * few instances
	 * 
	 * @param graph
	 * @param sub
	 * @param occurrences
	 * @param multiEdges An empty frequencymodel receiving how often certain edges in 
	 * the template graph should be repeated (one occurrence in the fm no repeats). 
	 * For performance reasons, the actual indices refer to the old graph, not 
	 * the template graph. 
	 * @param rewiring An empty list, receiving the sequence of rewiring integers.
	 * @return
	 */
	public static Pair<Long, Long> subbedERInstances(
			UGraph<?> graph, UGraph<?> sub, List<List<Integer>> occurrences,
			FrequencyModel<Pair<Integer, Integer>> multiEdges,
			List<List<Integer>> rewiring)
	{
		if(occurrences.isEmpty())
			return p((long)graph.size(), graph.numLinks());
		
		long subbedSize = graph.size() - occurrences.size() * (occurrences.get(0).size() - 1);
		long subbedNumLinks = graph.numLinks() - sub.numLinks() * occurrences.size();
		// - we still need to remove multiple links from subbedNumLinks
	
		// * Which nodes have been mapped to which instance node 
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		// * Links that have been rewired. We store the old link (ie. the left 
		//   and right indices in the original graph). 
		Set<Pair<Integer, Integer>> rewLinks = 
				new LinkedHashSet<Pair<Integer, Integer>>();
		
		for(List<Integer> occurrence : occurrences)
		{	
			List<Integer> rw = new LinkedList<Integer>(); 
			for(int index : occurrence)
				map.put(index, occurrence.get(0));
			
			// * Remove all links linking into an occurrence
			for(int i : series(occurrence.size()))
			{
				int index = occurrence.get(i);
				for(UNode<?> node : graph.get(index).neighbors())
					if(! occurrence.contains(node.index()))
					{		
						rewLinks.add(ordered(index, node.index()));
						rw.add(i);						
					}
			}
				
			rewiring.add(rw);
		}
		
		subbedNumLinks -= rewLinks.size();
				
		// * convert the rewritten links to new indices, and build a 
		//  frequencymodel 
		for(Pair<Integer, Integer> link : rewLinks)
		{
			int f = link.first(), s = link.second();
			int a = map.containsKey(f) ? map.get(f) : f;
			int b = map.containsKey(s) ? map.get(s) : s;
			
			multiEdges.add(ordered(a, b));
		}
		
		// * Add each rewritten link _once_
		for(Pair<Integer, Integer> link : multiEdges.tokens())
			subbedNumLinks ++;

		return Pair.p(subbedSize, subbedNumLinks);
	}
	
	/**
	 * Computes the size and number of links in the template graph, by looping 
	 * only over the instances. This method should be fast for large graphs with
	 * few instances
	 * 
	 * @param graph
	 * @param sub
	 * @param occurrences
	 * @param multiEdges An empty frequencymodel receiving how often certain edges in 
	 * the template graph should be repeated (one occurrence in the fm no repeats). 
	 * For performance reasons, the actual indices refer to the old graph, not 
	 * the template graph. 
	 * @param rewiring An empty list, receiving the sequence of rewiring integers.
	 * @return
	 */
	public static Pair<Long, Long> subbedERInstances(
			DGraph<?> graph, DGraph<?> sub, List<List<Integer>> occurrences,
			FrequencyModel<Pair<Integer, Integer>> multiEdges,
			List<List<Integer>> rewiring)
	{
		if(occurrences.isEmpty())
			return p((long)graph.size(), graph.numLinks());
		
		long subbedSize = graph.size() - occurrences.size() * (occurrences.get(0).size() - 1);
		long subbedNumLinks = graph.numLinks() - sub.numLinks() * occurrences.size();
		// - we still need to remove multiple links from subbedNumLinks
	
		// * Which nodes have been mapped to which instance node 
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		// * Links that have been rewired. We store the old link (ie. the left 
		//   and right indices in the original graph). 
		Set<Pair<Integer, Integer>> rewLinks = 
				new LinkedHashSet<Pair<Integer, Integer>>();
		
		for(List<Integer> occurrence : occurrences)
		{	 

			List<Integer> rw = new LinkedList<Integer>(); 
			for(int index : occurrence)
				map.put(index, occurrence.get(0));
			
			// * Remove all links linking into an occurrence
			for(int i : series(occurrence.size()))
			{
				int index = occurrence.get(i);
				for(DNode<?> node : graph.get(index).out())
					if(! occurrence.contains(node.index()))
					{		
						rewLinks.add(p(index, node.index()));
						
						rw.add(i);						
					}
				
				for(DNode<?> node : graph.get(index).in())
					if(! occurrence.contains(node.index()))
					{		
						rewLinks.add(p(node.index(), index));
						
						rw.add(i);						
					}
			}
				
			rewiring.add(rw);

		}
		
		subbedNumLinks -= rewLinks.size();
				
		// * convert the rewritten links to new indices, and build a 
		//  frequencymodel 
		for(Pair<Integer, Integer> link : rewLinks)
		{
			int f = link.first(), s = link.second();
			int a = map.containsKey(f) ? map.get(f) : f;
			int b = map.containsKey(s) ? map.get(s) : s;
						
			multiEdges.add(p(a, b));
		}
				
		// * Add each rewritten link _once_
		for(Pair<Integer, Integer> link : multiEdges.tokens())
			subbedNumLinks ++;

		return Pair.p(subbedSize, subbedNumLinks);
	}
	
	/**
	 * Returns a pair of ordered integer (ie. the smallest integer always takes 
	 * the first position).
	 *  
	 * @param i1
	 * @param i2
	 * @return
	 */
	private static Pair<Integer, Integer> ordered(int i1, int i2) 
	{
		if(i1 <= i2)
			return new Pair<Integer, Integer>(i1, i2);
		return new Pair<Integer, Integer>(i2, i1);
	}
}
