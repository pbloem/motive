package nl.peterbloem.motive;

import static nl.peterbloem.kit.Series.series;
import static org.nodes.motifs.MotifCompressor.exDegree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.nodes.DGraph;
import org.nodes.DLink;
import org.nodes.DNode;
import org.nodes.DTGraph;
import org.nodes.DTLink;
import org.nodes.DTNode;
import org.nodes.Graph;
import org.nodes.Graphs;
import org.nodes.MapDTGraph;
import org.nodes.Node;
import org.nodes.Subgraph;
import org.nodes.algorithms.Nauty;
import org.nodes.compression.EdgeListCompressor;
import org.nodes.compression.Functions;
import org.nodes.compression.NeighborListCompressor;
import org.nodes.compression.Functions.NaturalComparator;
import org.nodes.motifs.MotifCompressor;
import org.nodes.motifs.MotifVarTags;
import org.nodes.random.SimpleSubgraphGenerator;
import org.nodes.random.SubgraphGenerator;
import org.nodes.util.Compressor;

import au.com.bytecode.opencsv.CSVWriter;
import nl.peterbloem.kit.AbstractGenerator;
import nl.peterbloem.kit.BitString;
import nl.peterbloem.kit.FrequencyModel;
import nl.peterbloem.kit.Generator;
import nl.peterbloem.kit.Generators;
import nl.peterbloem.kit.Global;
import nl.peterbloem.kit.Order;
import nl.peterbloem.kit.Series;

/**
 * Extracts motifs from a DGraph<String> by sampling.
 * 
 * This matches including node labels. It returns the frequent labeled subgraphs. 
 * If labels are to be ignored, the graph should be blanked beforehand (see 
 * Graphs.blank()).
 * 
 * @author Peter
 *
 */
public class DPlainMotifExtractor<L extends Comparable<L>>
{
	public static final boolean SORT_BY_SCORE = false;
	
	private DGraph<L> data;
	private int samples;

	private Functions.NaturalComparator<L> comparator;

	private List<DGraph<L>> tokens;

	private Generator<Integer> intGen;
	
	private MotifVarTags mvTop = null;
	
	private FrequencyModel<DGraph<L>> fm;
	private Map<DGraph<L>, List<List<Integer>>> occurrences;
	
	private int minFreq;
	
	public DPlainMotifExtractor(
			DGraph<L> data,
			int numSamples,
			int minSize,
			int maxSize,
			int minFreq)
	{
		this.data = data;
		this.samples = numSamples;
		
		comparator = new Functions.NaturalComparator<L>();
		intGen = Generators.uniform(minSize, maxSize + 1);
		
		this.minFreq = minFreq; 
		
		run();
	}
	
	/**
	 * Extracts motifs of a given size.
	 * 
	 * @param data
	 * @param numSamples
	 * @param size
	 */
	public DPlainMotifExtractor(
			DGraph<L> data,
			int numSamples,
			int size, 
			int minFreq)
	{
		this.data = data;
		this.samples = numSamples;
		
		comparator = new Functions.NaturalComparator<L>();
		intGen = Generators.uniform(size, size+1);
		
		this.minFreq = minFreq; 
		
		run();
	}

	private void run()
	{
		if(samples > 10000) Global.log().info("Sampling motifs");		
		fm = new FrequencyModel<DGraph<L>>();

		// * The (overlapping) instances
		occurrences = new LinkedHashMap<DGraph<L>, List<List<Integer>>>();

		SimpleSubgraphGenerator gen = 
			new SimpleSubgraphGenerator(data, intGen);
		
//		SubgraphGenerator<L> gen = new SubgraphGenerator<L>(data, intGen);

		if(samples > 10000) Global.log().info("Start sampling.");
		int interval = samples / 20;
		
		for (int i : Series.series(samples))
		{
			if (i % interval == 0 && samples > 10000)
				Global.log().info("Samples finished: " + i + " (" + ((i*100)/samples) + "%)");

			List<Integer> indices = gen.generate();
			
			DGraph<L> sub = Subgraph.dSubgraphIndices(data, indices);

			// * Reorder nodes to canonical ordering
			Order canonical = Nauty.order(sub, comparator);
			sub = Graphs.reorder(sub, canonical);
			
			List<Integer> occurrence = canonical.apply(indices); 
			
			fm.add(sub);
			
			// * record the occurrence
			if (!occurrences.containsKey(sub))
				occurrences.put(sub, new ArrayList<List<Integer>>());

			occurrences.get(sub).add(occurrence);
		}
		
		if(samples > 10000) Global.log().info("Removing overlapping occurrences.");
		// * Remove overlapping occurrences 
		//   (keep the ones with the lowest exdegrees)
		FrequencyModel<DGraph<L>> newFm = 
				new FrequencyModel<DGraph<L>>();
		Map<DGraph<L>, List<List<Integer>>> newOccurrences = 
				new LinkedHashMap<DGraph<L>, List<List<Integer>>>();
		
		for(DGraph<L> sub : fm.tokens())
			if(fm.frequency(sub) >= minFreq)
			{
				int nl = (int)sub.numLinks();
				
				// * A map from nodes to occurrences containing them
				Map<DNode<L>, List<Occurrence>> map = 
					new LinkedHashMap<DNode<L>, List<Occurrence>>();
				
				// * A list of all occurrences, sorted by exDegree
				LinkedList<Occurrence> list = new LinkedList<Occurrence>();
				
				// - fill the map and list
				for(List<Integer> occurrence : occurrences.get(sub))
				{
					Occurrence occ = new Occurrence(occurrence);
					
					list.add(occ);
					
					for(int index : occurrence)
					{
						DNode<L> node = data.get(index);
						
						if(! map.containsKey(node))
							map.put(node, new ArrayList<Occurrence>());	 
						map.get(node).add(occ);
					}
				}
				
				Collections.sort(list);
				
				while(!list.isEmpty())
				{
					// * Find the first living occurrence, remove any dead ones
					Occurrence head = list.poll();
					
					if(head.alive()) // * register it as a viable occurrence
					{
						newFm.add(sub, SORT_BY_SCORE ? nl/(double)head.exDegree() : 1.0);
						
						if(! newOccurrences.containsKey(sub))
							newOccurrences.put(sub, new ArrayList<List<Integer>>());
						newOccurrences.get(sub).add(head.indices());
						
						// - now kill any occurrence that shares a node with this one.
						for(int nodeIndex : head.indices())
							for(Occurrence occ : map.get(data.get(nodeIndex)))
								occ.kill();
					} 
				}
			}	
		
		fm = newFm;
		occurrences = newOccurrences;
		
		if(samples > 10000) Global.log().info("Finished sampling motifs and removing overlaps.");

		tokens = new ArrayList<DGraph<L>>(fm.sorted());
	}
	
	public List<DGraph<L>> subgraphs()
	{
		return tokens;
	}
	
	public List<List<Integer>> occurrences(DGraph<L> subgraph)
	{
		return occurrences.get(subgraph);
	}

	private class Occurrence implements Comparable<Occurrence> {
		
		private boolean alive = true;
		private List<Integer> indices;
		private int id;
		
		private int exDegree;
	
		public Occurrence(List<Integer> indices)
		{
			this.indices = indices;
			
			exDegree = MotifCompressor.exDegree(data, indices);
		}
		
		public boolean alive()
		{
			return alive;
		}
		
		public void kill()
		{
			alive = false;
		}
		
		public List<Integer> indices()
		{
			return indices;
		}
		
		public int exDegree()
		{
			return exDegree;
		}

		@Override
		public int compareTo(Occurrence other)
		{
			return Integer.compare(exDegree, other.exDegree);
		}		
	}

	public double frequency(DGraph<L> sub)
	{
		return fm.frequency(sub);
	}
}
