package nl.peterbloem.motive.exec;

import static java.lang.Math.max;
import static nl.peterbloem.kit.Functions.log2;
import static nl.peterbloem.kit.Series.series;
import static org.apache.commons.math3.util.ArithmeticUtils.binomialCoefficientLog;
import static org.nodes.models.USequenceEstimator.CIMethod;
import static org.nodes.models.USequenceEstimator.CIType;
import static org.nodes.motifs.MotifCompressor.exDegree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

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
import static org.nodes.models.DSequenceEstimator.D;
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

public class Konect
{
	
	public static final double THRESHOLD = - log2(0.01);

	public String wgetprefix = "";
	public String tarprefix = "";
	/**
	 * Maximum amount of rewritten links.
	 */
	public int id = 0;

	/**
	 * Maximum amount of rewritten links.
	 */
	public boolean undirected = true;

	/**
	 * Maximum amount of rewritten links.
	 */
	public int maxRW = -1;

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
	public Graph<String> data = null;
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
	 * Depth to which to search for which instances to discard (more is better
	 * but slower, -1 is maximal depth always).
	 */
	public int searchDepth = -1;

	/**
	 * Whether to reset the DM model at every motif instance.
	 */
	private boolean resets = true;
	
	private double sampleTime = -1.0;

	public void main() throws IOException
	{
		Global.randomSeed();

		loadData();
		
		int nodes = data.size();
		long links = data.numLinks();
		
		int sig;
		long t0 = System.nanoTime();
		if(undirected)
			sig = undirected();
		else
			sig = directed();
		
		double elapsed = (System.nanoTime() - t0) / 1.0e9;
		
		String rand = String.format("%d05", Global.random().nextInt(10000));
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File("output."+rand+".csv")));
		writer.write(id + ", " + undirected + ", " + nodes + ", " + links + ", " + sig + ", " + sampleTime + ", " + elapsed + "\n");		
		writer.close();
	}

	public void loadData()
	{
		String url = "", filename = "";
		if (undirected)
		{
			switch (id)
			{
			case 0:
				dataName = "arenas-email";
				url = "http://konect.uni-koblenz.de/downloads/tsv/arenas-email.tar.bz2";
				filename = "arenas-email/out.arenas-email";
				break;
			case 1:
				dataName = "euroroad";
				url = "http://konect.uni-koblenz.de/downloads/tsv/subelj_euroroad.tar.bz2";
				filename = "subelj_euroroad/out.subelj_euroroad_euroroad";
				break;
			case 2:
				dataName = "chicago";
				url = "http://konect.uni-koblenz.de/downloads/tsv/tntp-ChicagoRegional.tar.bz2";
				filename = "tntp-ChicagoRegional/out.tntp-ChicagoRegional";
				break;
			case 3:
				dataName = "hamsterster-friendships";
				url = "http://konect.uni-koblenz.de/downloads/tsv/petster-friendships-hamster.tar.bz2";
				filename = "petster-friendships-hamster/out.petster-friendships-hamster-uniq";
				break;
			case 4:
				dataName = "hamsterster-full";
				url = "http://konect.uni-koblenz.de/downloads/tsv/petster-hamster.tar.bz2";
				filename = "petster-hamster/out.petster-hamster";
				break;
			case 5:
				dataName = "facebook-nips";
				url = "http://konect.uni-koblenz.de/downloads/tsv/ego-facebook.tar.bz2";
				filename = "ego-facebook/out.ego-facebook";
				break;
			case 6:
				dataName = "us-power-grid";
				url = "http://konect.uni-koblenz.de/downloads/tsv/opsahl-powergrid.tar.bz2";
				filename = "opsahl-powergrid/out.opsahl-powergrid";
				break;
			case 7:
				dataName = "arxiv-astro-ph";
				url = "http://konect.uni-koblenz.de/downloads/tsv/ca-AstroPh.tar.bz2";
				filename = "ca-AstroPh/out.ca-AstroPh";
				break;
			case 8:
				dataName = "brightkite";
				url = "http://konect.uni-koblenz.de/downloads/tsv/loc-brightkite_edges.tar.bz2";
				filename = "loc-brightkite_edges/out.loc-brightkite_edges";
				break;
			case 9:
				dataName = "livemocha";
				url = "http://konect.uni-koblenz.de/downloads/tsv/livemocha.tar.bz2";
				filename = "livemocha/out.livemocha";
				break;
			case 10:
				dataName = "flickr";
				url = "http://konect.uni-koblenz.de/downloads/tsv/flickrEdges.tar.bz2";
				filename = "flickrEdges/out.flickrEdges";
				break;
			case 11:
				dataName = "wordnet";
				url = "http://konect.uni-koblenz.de/downloads/tsv/wordnet-words.tar.bz2";
				filename = "wordnet-words/out.wordnet-words";
				break;
			case 12:
				dataName = "douban";
				url = "http://konect.uni-koblenz.de/downloads/tsv/douban.tar.bz2";
				filename = "douban/out.douban";
				break;
			case 13:
				dataName = "gowalla";
				url = "http://konect.uni-koblenz.de/downloads/tsv/loc-gowalla_edges.tar.bz2";
				filename = "loc-gowalla_edges/out.loc-gowalla_edges";
				break;
			case 14:
				dataName = "dblp";
				url = "http://konect.uni-koblenz.de/downloads/tsv/com-dblp.tar.bz2";
				filename = "com-dblp/out.com-dblp";
				break;
			case 15:
				dataName = "amazon-mds";
				url = "http://konect.uni-koblenz.de/downloads/tsv/com-amazon.tar.bz2";
				filename = "com-amazon/out.com-amazon";
				break;
			case 16:
				dataName = "pennsylvania";
				url = "http://konect.uni-koblenz.de/downloads/tsv/roadNet-PA.tar.bz2";
				filename = "roadNet-PA/out.roadNet-PA";
				break;
			case 17:
				dataName = "youtube-friendship";
				url = "http://konect.uni-koblenz.de/downloads/tsv/com-youtube.tar.bz2";
				filename = "com-youtube/out.com-youtube";
				break;
			}
		} else
		{
			switch (id)
			{
			case 0:
				dataName = "human-protein-figeys";
				url = "http://konect.uni-koblenz.de/downloads/tsv/maayan-figeys.tar.bz2";
				filename = "maayan-figeys/out.maayan-figeys";
				break;
			case 1:
				dataName = "openflights";
				url = "http://konect.uni-koblenz.de/downloads/tsv/opsahl-openflights.tar.bz2";
				filename = "opsahl-openflights/out.opsahl-openflights";
				break;
			case 2:
				dataName = "twitter-lists";
				url = "http://konect.uni-koblenz.de/downloads/tsv/ego-twitter.tar.bz2";
				filename = "ego-twitter/out.ego-twitter";
				break;
			case 3:
				dataName = "google-plus";
				url = "http://konect.uni-koblenz.de/downloads/tsv/ego-gplus.tar.bz2";
				filename = "ego-gplus/out.ego-gplus";
				break;
			case 4:
				dataName = "gnutella";
				url = "http://konect.uni-koblenz.de/downloads/tsv/p2p-Gnutella31.tar.bz2";
				filename = "p2p-Gnutella31/out.p2p-Gnutella31";
				break;
			case 5:
				dataName = "epinions";
				url = "http://konect.uni-koblenz.de/downloads/tsv/soc-Epinions1.tar.bz2";
				filename = "soc-Epinions1/out.soc-Epinions1";
				break;
			case 6:
				dataName = "stanford";
				url = "http://konect.uni-koblenz.de/downloads/tsv/web-Stanford.tar.bz2";
				filename = "web-Stanford/out.web-Stanford";
				break;
			case 7:
				dataName = "amazon-tweb";
				url = "http://konect.uni-koblenz.de/downloads/tsv/amazon0601.tar.bz2";
				filename = "amazon0601/out.amazon0601";
				break;
			case 8:
				dataName = "twitter-icwsm";
				url = "http://konect.uni-koblenz.de/downloads/tsv/munmun_twitter_social.tar.bz2";
				filename = "munmun_twitter_social/out.munmun_twitter_social";
				break;
			case 9:
				dataName = "berkeley-stanford";
				url = "http://konect.uni-koblenz.de/downloads/tsv/web-BerkStan.tar.bz2";
				filename = "web-BerkStan/out.web-BerkStan";
				break;
			case 10:
				dataName = "google";
				url = "http://konect.uni-koblenz.de/downloads/tsv/web-Google.tar.bz2";
				filename = "web-Google/out.web-Google";
				break;
			case 11:
				dataName = "youtube-links";
				url = "http://konect.uni-koblenz.de/downloads/tsv/youtube-links.tar.bz2";
				filename = "youtube-links/out.youtube-links";
				break;
			}
		}

		try
		{
			Process p;
			p = Runtime.getRuntime().exec(wgetprefix + "wget " + url);
			print(p);

			p = Runtime.getRuntime().exec(new String[]
			{ "/bin/bash", "-c", tarprefix + "tar xvfj *.tar.bz2" });
			print(p);

			if (undirected)
				data = Data.edgeListUndirectedUnlabeled(new File("./" + filename), true);
			else
				data = Data.edgeListDirectedUnlabeledSimple(new File("./" + filename));

		} catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}

	private void print(Process p)
		throws IOException
	{
		String s= null;
		BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));

		BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));

		// read the output from the command
		while ((s = stdInput.readLine()) != null)
			System.out.println(s);

		// read any errors from the attempted command
		while ((s = stdError.readLine()) != null)
			System.out.println(s);
	}
	
	public int directed()
	{
		
		int significant = 0;
		
		DGraph<String> data = (DGraph<String>) this.data;

		List<D> degrees = DSequenceEstimator.sequence((DGraph<String>) data);

		long t0 = System.nanoTime();
		DPlainMotifExtractor<String> ex = new DPlainMotifExtractor<String>(data, motifSamples,
				motifMinSize, motifMaxSize, minFreq);
		this.sampleTime = (System.nanoTime() - t0)/1.0e9;
		

		List<? extends DGraph<String>> subsAll = new ArrayList<DGraph<String>>(ex.subgraphs());
		List<Double> frequenciesAll = new ArrayList<Double>(subsAll.size());

		for (DGraph<String> sub : subsAll)
			frequenciesAll.add(ex.frequency(sub));

		List<List<List<Integer>>> occurrences = new ArrayList<List<List<Integer>>>(subsAll.size());
		for (DGraph<String> sub : subsAll)
			occurrences.add(ex.occurrences(sub));

		// - select the top motifs by frequency
		List<? extends DGraph<String>> subs;
		List<Double> frequencies;
		if (subsAll.size() > maxMotifs)
		{
			subs = new ArrayList<DGraph<String>>(subsAll.subList(0, maxMotifs));
			frequencies = new ArrayList<Double>(frequenciesAll.subList(0, maxMotifs));
		} else
		{
			subs = subsAll;
			frequencies = frequenciesAll;
		}

		double baselineEL = EdgeListModel.directed(degrees, Prior.ML);

		for (final int i : series(subs.size()))
		{
			DGraph<String> sub = subs.get(i);
			List<List<Integer>> occs = occurrences.get(i);

			double max = Double.NEGATIVE_INFINITY;

			double sizeEL = MotifSearchModel.sizeELInst(data, degrees, sub, occs, resets, searchDepth);
			double factorEL = baselineEL - sizeEL;
			
			if(factorEL > THRESHOLD)
				significant ++;
		}
		
		return significant;
	}
	
	public int undirected()
	{
		int significant = 0;
		
		UGraph<String> data = (UGraph<String>) this.data;

		List<Integer> degrees = Graphs.degrees(data);

		long t0 = System.nanoTime();
		UPlainMotifExtractor<String> ex = new UPlainMotifExtractor<String>(data, motifSamples,
				motifMinSize, motifMaxSize, minFreq);
		this.sampleTime = (System.nanoTime() - t0)/1.0e9;

		List<? extends UGraph<String>> subsAll = new ArrayList<UGraph<String>>(ex.subgraphs());
		List<Double> frequenciesAll = new ArrayList<Double>(subsAll.size());

		for (UGraph<String> sub : subsAll)
			frequenciesAll.add(ex.frequency(sub));

		final List<List<List<Integer>>> occurrences = new ArrayList<List<List<Integer>>>(subsAll.size());
		for (UGraph<String> sub : subsAll)
			occurrences.add(ex.occurrences(sub));

		// - select the top motifs by frequency
		List<? extends UGraph<String>> subs;
		List<Double> frequencies;
		if (subsAll.size() > maxMotifs)
		{
			subs = new ArrayList<UGraph<String>>(subsAll.subList(0, maxMotifs));
			frequencies = new ArrayList<Double>(frequenciesAll.subList(0, maxMotifs));
		} else
		{
			subs = subsAll;
			frequencies = frequenciesAll;
		}

		double baselineEL = EdgeListModel.undirected(degrees, Prior.ML);

		for (final int i : series(subs.size()))
		{
			UGraph<String> sub = subs.get(i);
			List<List<Integer>> occs = occurrences.get(i);

			double max = Double.NEGATIVE_INFINITY;

			double sizeEL = MotifSearchModel.sizeELInst(data, degrees, sub, occs, resets, searchDepth);
			double factorEL = baselineEL - sizeEL;
			
			if(factorEL > THRESHOLD)
				significant ++;
		}
		
		return significant;
	}
}
