package nl.peterbloem.motive.exec;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.nodes.compression.Functions;

import nl.peterbloem.kit.Global;

public class Run {

	
	@Option(
		name="--type",
		usage="Selects the type of experiment, one of: synth (synthetic graph experiment), full (motif extraction with all null models), fast (skip the DS model)}.")
	private static String type = "fast";
	
	private static int minSize = 3;
	
	private static int maxSize = 6;
	
	@Option(name="--file", usage="Input file: a graph in edge-list encoding (2 tab separated integers per line)")
	private static File graph;
	
	@Option(name="--directed", usage="Whether the input is directed or not.")
	private static boolean directed;
	
	
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
    	
    	if ("synth".equals(type.toLowerCase()))
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
    		
    		
    	} else if ("full".equals(type.toLowerCase()))
    	{
    		Global.log().info("Experiment type: full");
		
    		
    	} else 
    	{
    		Global.log().severe("Experiment type " + type + " not recognized. Exiting.");
    		System.exit(-1);
    	}
    	
    	

	}

}
