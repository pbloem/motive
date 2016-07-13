# motive
A proof-of-concept library for motif analysis using MDL techniques. It contains the methods described in the paper _Compression as a Fast Measure of Network Motif Relevance_.

## Installation and use

Each release on github comes with a compiled JAR file which you can run as a command-line program. Download [the latest one here](https://github.com/pbloem/motive/releases). See the section _examples_ below, for how to use it.

If you would like to call the code directly from your own program, the simplest way is to include it as a Maven dependency through [jitpack](http://jitpack.io/pbloem/motive). Just include the following repository in your pom:

```xml
    <repositories>
        <repository>
            <id>jitpack.io</id>
            <url>https://jitpack.io</url>
        </repository>
    </repositories>
```

and the following dependency:

```xml
	<dependency>
	    <groupId>com.github.pbloem</groupId>
	    <artifactId>motive</artifactId>
	    <version>v0.1.0</version>
	</dependency>
```
Check the jitpack link above for linking from gradle/sbt/leiningen projects, and to see what the latest release is.

Have a look at the classes Compare and CompareLarge for hints on how to set up a motif experiment from within java code. For comman line usage, see the next section.

## Usage examples

Display usage information:
```bash
java -jar motive.jar --help 
```

Run the "synthetic" experiment: create random graphs with injected motifs.

```bash
java -jar motive.jar --type synth --synth.repeats 3 --synth.n 50 --synth.m 600 --synth.instances 0,5,10
```

Run the "fast" experiment (ER model and EdgeList model) on a particular dataset, finding motifs of size up to (and including) 10:

```bash
java -jar motive.jar --type fast --file /Users/Peter/Documents/Datasets/graphs/cit/simple.txt --minsize 3 --maxsize 10 --samples 1000000 --maxmotifs 30 
```

The "full" experiment includes the precise DS model as well. This is a bit slower.
```bash
java -jar motive.jar --type full --file /Users/Per/Documents/Datasets/graphs/cit/simple.txt --minsize 3 --maxsize 5 --samples 100000 --maxmotifs 30

```  
## Plotting

The plots in the paper were produced using python scripts. Motive will copy these into the output directory and attempt to run them. This will fail if the correct dependencies are not installed. Here's a short recipe for getting the scripts to run:

 * Install [Anaconda](https://www.continuum.io/downloads) (the python 3 version)
 * run the command ```conda install matplotlib'''
 * Go to the output directory and run ```python plot.synthetic.py''' (or whichever python file was copied to the output directory).
 
If this doesn't work, [make a ticket](https://github.com/pbloem/motive/issues) on github describing what went wrong, and we'll try to help you out. 

## Miscelaneous notes

* The Nauty implementation, used for graph canonization is a re-implementation. It is not complete, and the real Nauty is likely much faster. If many samples of large motifs (eg. 12 nodes) are necessary, this will become the bottleneck. Up to 10 nodes, however, you should be able to take around a 25k 10-node samples per minute, even on commodity hardware. 
* If you get a memory error, increase your heap space by adding ```-Xmx3000m''' before "-jar". This sets the heap space to 3 gigabytes. Don't set the heap space to more memory than your machine has available (and leave at least  500m to 1000m for the system). If you still get memory errors, reduce the maximum motif size, the number of samples, or see if you can work with a smaller dataset (or move to a machine with more memory, of course).
 