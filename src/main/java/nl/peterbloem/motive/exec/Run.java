package nl.peterbloem.motive.exec;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.nodes.DGraph;
import org.nodes.DiskDGraph;
import org.nodes.Graph;
import org.nodes.compression.Functions;
import org.nodes.data.Data;
import org.nodes.data.GML;
import org.nodes.data.RDF;

import nl.peterbloem.kit.Global;

public class Run 
{
	
	@Option(
		name="--type",
		usage="Selects the type of experiment, one of: synth (synthetic graph experiment), full (motif extraction with all null models), fast (skip the DS model), preload (load a large graph into a db file), class (classification experiment).")
	private static String type = "fast";
	
	@Option(
			name="--preload.out",
			usage="Output file.")
	private static File outFile = new File("./graph.db");
	
	@Option(
			name="--samples",
			usage="Number of samples to take.")
	private static int samples = 1000000;
	
	@Option(
			name="--minsize",
			usage="Minimum motif size in nodes (inclusive).")
	private static int minSize = 3;
	@Option(
			name="--maxsize",
			usage="Maximum motif size in nodes (inclusive).")
	private static int maxSize = 6;
	
	@Option(
			name="--maxmotifs",
			usage="Maximum number of motifs to test.")
	private static int maxMotifs = 100;
	
	@Option(name="--file", usage="Input file: a graph in edge-list encoding (2 tab separated integers per line). Multiple edges and self-loops are ignored. If type is class, this should be an RDF file.")
	private static File file;
	
	@Option(name="--class.table", usage="TSV file containing the classification experiment.")
	private static File classTSV;
	
	@Option(name="--filetype", usage="Filetype: edgelist, gml or graph (ie. a file created with \"--type preload\")")
	private static String filetype = "edgelist";
	
	@Option(name="--undirected", usage="If the input should be interpeted as undirected (only for edgelist files).")
	private static boolean undirected = false;
	
	@Option(name="--threads", usage="Number of threads to run simultaneaously. Default is the number of cores available. In full mode, there will always be at least 2 concurrent threads, even if this value is 1.")
	private static int threads = Global.numThreads();
	
	@Option(
			name="--fast.graphloop",
			usage="Loop over the graph instead of the instances when computing the score. A little faster when there are many instances, but a lot slower when there are few.")
	private static boolean graphLoop = false;
	
	@Option(
			name="--fast.disk",
			usage="Use the disk to store the graph.  Slower, but uses very little memory. Supports graphs up to billions of links (disk space permitting).")
	private static boolean useDisk = false;
	
	@Option(
			name="--full.depth",
			usage="The search depth for the DS model.")
	private static int dsDepth = 3;
	
	@Option(
			name="--full.mix",
			usage="What proportion of available cores to use for motif computation (with the rest used for sampling). Changing this parameter won't affect the end result, but it might lead to better utilitzation of the available cores. If 1.0, the motif scores are computed one by one, sequentially, and the free cores are used to sample for the DS model. If 0.0, the motif scores are computed in parallel, and sampling is done single-threaded.")
	private static double mix = 0.4;
	
	@Option(name="--help", usage="Print usage information.", aliases={"-h"}, help=true)
	private static boolean help = false;
	
	@Option(
			name="--synth.repeats",
			usage="How often the synth experiment should be repeated.")
	private static int synthRepeats = 10;
	
	@Option(
			name="--synth.n",
			usage="Number of nodes in the sampled graph.")
	private static int synthN = 5000;
	
	@Option(
			name="--synth.m",
			usage="Number of links in the sampled graph.")
	private static int synthM = 10000;
	
	@Option(
			name="--synth.instances",
			usage="Number of motif instances to inject. Should be a list of comma-separated ints (no spaces).")
	private static String synthNumInstances = "0,10,100";
	
	@Option(
			name="--synth.maxdegree",
			usage="Maximum degree for an instance node.")
	private static int synthMaxDegree = 5;
	
	@Option(
			name="--class.prob",
			usage="Fanmod probability of expanding a search tree node (1.0 enumerates all subgraphs, lower means a smaller sample, and lower runtime)")
	private static double classProb = 0.5;
	
	@Option(
			name="--class.hubs",
			usage="Number of hubs to remove from the data (the more hubs removed, the smaller the instances become.")
	private static int classHubs = 0;
	
	@Option(
			name="--class.fanmodSamples",
			usage="Number of samples from the null model in the FANMOD experiment.")
	private static int classFMSamples = 1000;
	
	@Option(
			name="--class.motiveSamples",
			usage="Number of subgraphs to sample in the motive experiment.")
	private static int classMotiveSamples = 1000000;
	
	@Option(
			name="--class.depth",
			usage="Depth to which to extract the instances.")
	private static int classDepth = 2;
	
	@Option(
			name="--class.mixingTime",
			usage="Mixing time for the curveball sampling algorithm (ie. the number of steps taken in the markov chain for each sample).")
	private static int classMixingTime = 10000;
	
	@Option(
			name="--class.numInstances",
			usage="The number of instances to use (samples from the total available)")
	private static int classNumInstances = 100;
	
	@Option(
			name="--class.sizes",
			usage="The motif sizes to use as features")
	private static String classSizes = "3,4";
	
	/**
	 * Main executable function
	 * @param args
	 */
	public static void main(String[] args) 
	{	
		Run run = new Run();
		
		// * Parse the command-line arguments
    	CmdLineParser parser = new CmdLineParser(run);
    	try
		{
			parser.parseArgument(args);
		} catch (CmdLineException e)
		{
	    	System.err.println(e.getMessage());
	        System.err.println("java -jar motive.jar [options...]");
	        parser.printUsage(System.err);
	        
	        System.exit(1);	
	    }
    	
    	if(help)
    	{
	        parser.printUsage(System.out);
	        
	        System.exit(0);	
    	}
    	
		Global.setNumThreads(threads);
		Global.log().info("Using " + Global.numThreads() + " concurrent threads");
    	
    	if ("class".equals(type.toLowerCase()))
    	{
    	
    		ClassExperiment exp = new ClassExperiment();
    		
    		try {
				exp.graph = RDF.readSimple(file);
			} catch (IOException e) {
				throw new RuntimeException("Could not read RDF input file.", e);
			}
    		
    		try {
    			exp.map = ClassExperiment.tsv(classTSV);
			} catch (IOException e) {
				throw new RuntimeException("Could not read TSV classification file.", e);
			}	
    		
    		exp.prob = classProb;
    		exp.hubsToRemove = classHubs;
    		exp.samples = classFMSamples;
    		exp.motiveSamples = classMotiveSamples;
    		exp.instanceDepth = classDepth;
    		exp.mixingTime = classMixingTime;
    		exp.numInstances = classNumInstances;

    		exp.sizes = new ArrayList<Integer>();
    		try{
	    		for(String elem : classSizes.split(","))
	    		{
	    			exp.sizes.add(Integer.parseInt(elem));
	    		}
    		} catch(RuntimeException e)
    		{
    			throw new RuntimeException("Failed to parse sizes argument: " + classSizes + " (does it contain spaces, or non-integers?)." , e);
    		}
    		
    		try {
				exp.main();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
    		    		
    		
    	} else if ("preload".equals(type.toLowerCase()))
    	{
    		try 
    		{    			
    			File tmpDir = new File("./tmp/");
    			tmpDir.mkdir();
    			
    			DiskDGraph graph = DiskDGraph.fromFile(file, tmpDir, outFile);
    			graph.close();
    			
    			for(File file : tmpDir.listFiles())
    				file.delete();
    			tmpDir.delete();
    		} catch(IOException e)
    		{
    			throw new RuntimeException("Something went wrong reading or writing files.", e);
    		}
    		
    	} else if ("synth".equals(type.toLowerCase()))
    	{
    		Global.log().info("Experiment type: synthetic");
    		Synthetic synth = new Synthetic();
    		
    		// * Set the parameters
    		synth.runs = synthRepeats;
    		synth.n = synthN;
    		synth.m = synthM;
    		synth.maxDegree = synthMaxDegree;
    		
    		List<Integer> numInstances = new ArrayList<Integer>();
    		try{
	    		for(String elem : synthNumInstances.split(","))
	    		{
	    			numInstances.add(Integer.parseInt(elem));
	    		}
    		} catch(RuntimeException e)
    		{
    			throw new RuntimeException("Failed to parse num instances argument: " + synthNumInstances + " (does it contain spaces, or non-integers?)." , e);
    		}
    		
    		synth.numsInstances = numInstances;
    		
    		// * Run the experiment
    		Global.log().info("Starting experiment.");
    		Functions.tic();
    		try {
    			synth.main();
    		} catch(IOException e)
    		{
    			throw new RuntimeException("Encountered a problem when writing the results to disk.", e);
    		}
    		
    		Global.log().info("Experiment finished. Time taken: "+(Functions.toc())+" seconds.");
    	} else if ("fast".equals(type.toLowerCase()))
    	{
    		Global.log().info("Experiment type: fast");
    		
    		if(undirected)
    			throw new IllegalArgumentException("fast experiment is currently only implemented for directed graphs (please make a ticket on the github page if you need this feature for undirected graphs).");
    		
    		DGraph<String> data;
    		try {
	    		if("edgelist".equals(filetype.toLowerCase().trim()))
	    		{
	    			data = useDisk ? DiskDGraph.fromFile(file, new File("./tmp/")) : Data.edgeListDirectedUnlabeledSimple(file);
	    		} else if ("gml".equals(filetype.toLowerCase().trim()))
	    		{
	    			if(!useDisk)
	    				throw new IllegalArgumentException("GML is not supported with mode fast.disk");
	    			
	    			Graph<String> graph = GML.read(file);
	    			
	    			if(! (graph instanceof DGraph<?>))
	    				throw new IllegalArgumentException("Input file seems to describe an undirected graph. This is not (yet) supported for the fast mode.");
	    			
	    			data = (DGraph<String>) graph;
	    		} else if("db".equals(filetype.toLowerCase().trim())) 
	    		{
	    			data = DiskDGraph.fromDB(file);
    			} else {
	    			throw new IllegalArgumentException("File type ("+filetype+") not recognized.");
	    		}
	    			
    			
    		} catch (IOException e) {
				throw new IllegalArgumentException("There was a problem reading the input file ("+file+").", e);
			}
    		
    		CompareLarge large = new CompareLarge();
    		
    		large.dataName = file.getName();
    		large.data = data;
    		large.motifMinSize = minSize;
    		large.motifMaxSize = maxSize;
    		large.maxMotifs = maxMotifs;
    		large.motifSamples = samples;
    		large.graphLoop = graphLoop;
    		
       		Global.log().info("Starting experiment.");
    		Functions.tic();
    		try {
    			large.main();
    		} catch(IOException e)
    		{
    			throw new RuntimeException("Encountered a problem when writing the results to disk.", e);
    		}
    		
    		Global.log().info("Experiment finished. Time taken: "+(Functions.toc())+" seconds.");
    		
    		if(useDisk)
    		{
    			File dir = new File("./tmp/");
    			for(File file : dir.listFiles())
    				file.delete();
    			dir.delete();
    		}
    		
    	} else if ("full".equals(type.toLowerCase()))
    	{
    		Global.log().info("Experiment type: full");
		
    		Graph<String> data;
    		
    		
    		try {
    			if("edgelist".equals(filetype.trim().toLowerCase()))
    			{
	    			if(undirected)
	    				data = Data.edgeListUndirectedUnlabeled(file, true); // TODO: read as simple graph (like directed)
	    			else
	    				data = Data.edgeListDirectedUnlabeledSimple(file);
    			} else if ("edgelist".equals(filetype.trim().toLowerCase()))
    			{
    				data = GML.read(file);
    			} else {
	    			throw new IllegalArgumentException("File type ("+filetype+") not recognized.");
	    		}
			} catch (IOException e) {
				throw new IllegalArgumentException("There was a problem reading the input file ("+file+").");
			}
    		
    		Compare full = new Compare();
    		
    		full.dataName = file.getName();
    		full.data = data;
    		full.motifMinSize = minSize;
    		full.motifMaxSize = maxSize;
    		full.maxMotifs = maxMotifs;
    		full.motifSamples = samples;
    		full.betaSearchDepth = dsDepth;
    		full.mix = mix;
    		
       		Global.log().info("Starting experiment.");
    		Functions.tic();
    		try {
    			full.main();
    		} catch(IOException e)
    		{
    			throw new RuntimeException("Encountered a problem when writing the results to disk.", e);
    		}
    		
    		Global.log().info("Experiment finished. Time taken: "+(Functions.toc())+" seconds.");
    	} else 
    	{
    		Global.log().severe("Experiment type " + type + " not recognized. Exiting.");
	        System.err.println("java -jar motive.jar [options...]");
	        parser.printUsage(System.err);
    		System.exit(-1);
    	}
	}

}
