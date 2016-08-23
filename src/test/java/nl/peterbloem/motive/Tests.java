package nl.peterbloem.motive;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.junit.Test;
import org.nodes.DGraph;
import org.nodes.DNode;
import org.nodes.DiskDGraph;
import org.nodes.random.RandomGraphs;

/**
 * General tests, not specific to a class.
 * @author Peter
 *
 */
public class Tests {
	public static final File DIR = new File("./tmp/");
	
	@Test
	public void testSubgraphs()
		throws IOException
	{
		
		DGraph<String> graph = RandomGraphs.randomDirectedFast(100, 200);
		graph = DiskDGraph.copy(graph, DIR);
		
		DPlainMotifExtractor<String> ex = new DPlainMotifExtractor<String>(graph, 100000, 3, 6, 5);
		
		for(DGraph<String> sub : ex.subgraphs())
			for(DNode<String> node : sub.nodes())
			{
				assertTrue(node.neighbors().size() == node.degree());
				
				assertTrue(node.degree() > 0);
			}
	}

}
