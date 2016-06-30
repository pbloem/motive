# motive
A proof-of-concept library for motif analysis using MDL techniques. It contains the methods described in the paper _Compression as a Fast Measure of Network Motif Relevance_.

## Installation and use

Each release on github comes with a compiled JAR file which you can run as a command-line program. See the section _examples_ below.

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

## Examples

Run the "synthetic" experiment: create random graphs with injected motifs.

```bash
java -jar motive.jar --type synth --synth.repeats 3 --synth.n 50 --synth.m 600 --synth.instances 0,5,10
```

## Plotting

The plots in the paper were produced using python scripts. Motive will copy these into the output directory and attempt to run them. This will fail if the correct dependencies are not installed. Here's a short recipe for getting the scripts to run:

 * Install [Anaconda](https://www.continuum.io/downloads) (the python 3 version)
 * run the command ```conda install matplotlib'''
 * Go to the output directory and run ```python plot.synthetic.py''' (or whichever python file was copied to the output directory).
 
If this doesn't work, [make a ticket](https://github.com/pbloem/motive/issues) on github describing what went wrong, and we'll try to help you out. 
 