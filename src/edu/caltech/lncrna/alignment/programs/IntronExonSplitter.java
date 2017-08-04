package edu.caltech.lncrna.alignment.programs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Stream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.caltech.lncrna.bio.alignment.Aligned;
import edu.caltech.lncrna.bio.alignment.Alignment;
import edu.caltech.lncrna.bio.alignment.CoordinateSpace;
import edu.caltech.lncrna.bio.annotation.Annotated;
import edu.caltech.lncrna.bio.annotation.Annotation;
import edu.caltech.lncrna.bio.annotation.Annotation.AnnotationBuilder;
import edu.caltech.lncrna.bio.annotation.Strand;
import edu.caltech.lncrna.bio.datastructures.GenomeTree;
import edu.caltech.lncrna.bio.io.BamParser;
import edu.caltech.lncrna.bio.io.BamWriter;
import edu.caltech.lncrna.bio.io.BedParser;
import edu.caltech.lncrna.bio.io.SingleReadBamParser;

public final class IntronExonSplitter implements AutoCloseable {

    private final Path geneBedPath;
    private final Path bamPath;
    
    private final Path intronPath;
    private final Path exonPath;
    private final Path unclassifiedPath;
    
    private final BamParser<? extends Aligned<? extends Alignment>> bamParser;
    private final BamWriter intronWriter;
    private final BamWriter exonWriter;
    private final BamWriter unclassifiedWriter;
    
    private final static Path DEFAULT_INTRON_PATH = Paths.get("introns.bam");
    private final static Path DEFAULT_EXON_PATH = Paths.get("exons.bam");
    private final static Path DEFAULT_UNCLASSIFIED_PATH = Paths.get("unclassified.bam");
    private final static Path DEBUG_GENE_BODIES_PATH = Paths.get("gene_bodies.debug.bed");
    private final static Path DEBUG_INTRONS_PATH = Paths.get("introns.debug.bed");
    private final static int DEFAULT_EXON_PADDING = 0;
    
    private final CoordinateSpace header;
    
    private final GenomeTree<Annotated> introns = new GenomeTree<>();
    private final GenomeTree<Annotated> genes = new GenomeTree<>();
    private final GenomeTree<Annotated> paddedGenes = new GenomeTree<>();
    private final GenomeTree<Annotated> geneBodies = new GenomeTree<>();
    
    private final int exonPadding;
    private final boolean stranded;
    private final boolean forceSingle;
    private final boolean debug;
    private final static String VERSION = "1.0.0";
    
    private static final Logger LOGGER =
            Logger.getLogger(IntronExonSplitter.class.getName());
    
    public static void main(String[] args) throws IOException {

        CommandLine cmd = parseArgs(args);
        
        LOGGER.info("Starting process.");
        long startTime = System.currentTimeMillis();
        
        try (IntronExonSplitter splitter = new IntronExonSplitter(cmd)) {
            splitter
                .logParameters()
                .loadGenesAndIntrons()
                .getAlignments()
                .forEach(x -> splitter.writeAlignment(x));
        }
        
        long endTime = System.currentTimeMillis();
        LOGGER.info("Process complete. Took " + (endTime - startTime) + " ms.");
    }
    
    /**
     * Constructor.
     * <p>
     * Sets all variables as specified on the command line
     * @param cmd - the CommandLine object
     */
    private IntronExonSplitter(CommandLine cmd) {
        geneBedPath = Paths.get(cmd.getOptionValue("genes"));
        bamPath = Paths.get(cmd.getOptionValue("bam"));
        header = new CoordinateSpace(bamPath);
 
        intronPath = cmd.hasOption("introns")
                ? Paths.get(cmd.getOptionValue("introns"))
                : DEFAULT_INTRON_PATH;
        intronWriter = new BamWriter(intronPath, header);

        exonPath = cmd.hasOption("exons")
                ? Paths.get(cmd.getOptionValue("exons"))
                : DEFAULT_EXON_PATH;
        exonWriter = new BamWriter(exonPath, header);
        
        unclassifiedPath = cmd.hasOption("unclassified")
                ? Paths.get(cmd.getOptionValue("unclassified"))
                : DEFAULT_UNCLASSIFIED_PATH;
        unclassifiedWriter = new BamWriter(unclassifiedPath, header);

        debug = cmd.hasOption("debug");
        if (debug) {
            LOGGER.setLevel(Level.FINEST);
        }
        
        forceSingle = cmd.hasOption("single");
        
        bamParser = forceSingle
                ? new SingleReadBamParser(bamPath)
                : BamParser.newInstance(bamPath);
                
        exonPadding = cmd.hasOption("exonpadding")
                ? Integer.parseInt(cmd.getOptionValue("exonpadding"))
                : DEFAULT_EXON_PADDING;
                
        stranded = cmd.hasOption("stranded");
    }
  
    /**
     * Parses the command-line arguments.
     * <p>
     * If passed only "--help" or only "--version", prints a help menu or the
     * version number respectively, then exits.
     * <p>
     * Any parse errors cause the program to print a help menu, then exit with
     * a 1 status code.
     * @param args - the command-line arguments
     * @return a CommandLine object containing the argument information
     */
    private static CommandLine parseArgs(String[] args) {
        
        Option versionOption = Option.builder("v")
                .longOpt("version")
                .desc("print version and exit")
                .hasArg(false)
                .required(false)
                .build();
        
        Option helpOption = Option.builder("h")
                .longOpt("help")
                .desc("print usage information and exit")
                .hasArg(false)
                .required(false)
                .build();

        Option bamOption = Option.builder()
                .longOpt("bam")
                .desc("input BAM file")
                .argName("BAM")
                .hasArg()
                .required()
                .build();

        Option bedOption = Option.builder()
                .longOpt("genes")
                .desc("BED file of genes or transcripts")
                .argName("BED")
                .hasArg()
                .required()
                .build();
        
        Option exonOption = Option.builder()
                .longOpt("exons")
                .desc("output BAM file of reads contained entirely within " +
                        "exons; defaults to \"" +
                        DEFAULT_EXON_PATH.getFileName() + "\"")
                .argName("FILE")
                .hasArg()
                .required(false)
                .build();
        
        Option intronOption = Option.builder()
                .longOpt("introns")
                .desc("output BAM file of reads overlapping introns; " +
                        "defaults to \"" +
                        DEFAULT_INTRON_PATH.getFileName() + "\"")
                .argName("FILE")
                .hasArg()
                .required(false)
                .build();

        Option otherOption = Option.builder()
                .longOpt("unclassified")
                .desc("output BAM file of reads not overlapping any gene " + 
                        "body; defaults to \"" +
                        DEFAULT_UNCLASSIFIED_PATH.getFileName() + "\"")
                .argName("FILE")
                .hasArg()
                .required(false)
                .build();
        
        Option singleOption = Option.builder()
                .longOpt("single")
                .desc("force reading the BAM file as single-read alignments," +
                        " even if BAM file is paired-end.")
                .hasArg(false)
                .required(false)
                .build();
        
        Option debugOption = Option.builder()
                .longOpt("debug")
                .desc("run in debug mode")
                .hasArg(false)
                .required(false)
                .build();
        
        Option exonPaddingOption = Option.builder()
                .longOpt("exonpadding")
                .desc("pad exons by <NUM> bases")
                .argName("NUM")
                .hasArg()
                .required(false)
                .build();
        
        Option strandedOption = Option.builder()
                .longOpt("stranded")
                .desc("consider strandedness when calculating overlap")
                .hasArg(false)
                .required(false)
                .build();
        
        Options helpOptions = new Options().addOption(helpOption);
        Options versionOptions = new Options().addOption(versionOption);
        
        Options mainOptions = new Options()
                .addOption(bamOption)
                .addOption(bedOption)
                .addOption(intronOption)
                .addOption(exonOption)
                .addOption(otherOption)
                .addOption(debugOption)
                .addOption(singleOption)
                .addOption(exonPaddingOption)
                .addOption(strandedOption);

        Options allOptions = new Options();
        helpOptions.getOptions().forEach(allOptions::addOption);
        versionOptions.getOptions().forEach(allOptions::addOption);
        mainOptions.getOptions().forEach(allOptions::addOption);
        HelpFormatter formatter = new HelpFormatter();
                
        CommandLineParser parser = new DefaultParser();
        CommandLine rtrn = null;
        
        try {
            CommandLine cmds = parser.parse(helpOptions, args, true);
            if (cmds.getOptions().length == 1) {
                printHelp(formatter, allOptions);
                System.exit(0);
            }
            cmds = parser.parse(versionOptions, args, true);
            if (cmds.getOptions().length == 1) {
                System.out.println(VERSION);
                System.exit(0);
            }
            rtrn = parser.parse(mainOptions, args);
        } catch (ParseException e) {
            printHelp(formatter, allOptions);
            System.exit(1);
        }
        
        if (rtrn == null) {
            LOGGER.severe("An unknown error occurred while parsing command line arguments");
            System.exit(1);
        }
        
        return rtrn;
    }
    
    /**
     * Creates and prints a help menu to the console.
     * @param opts - the options to include in the help menu
     */
    private static void printHelp(HelpFormatter formatter, Options opts) {
        formatter.setDescPadding(0);
        String header = "\n";
        String footer = "\n";
        formatter.printHelp("java -jar IntronExonSplitter.jar", header,
                opts, footer, true);
    }
    
    /**
     * Log the given command-line parameters
     * @return this IntronExonSplitter for method chaining
     */
    private IntronExonSplitter logParameters() {
        LOGGER.info("BAM file of reads: " + bamPath.toString());
        LOGGER.info("BED file of gene annotations: " + geneBedPath.toString());
        LOGGER.info("Intron output file: " + intronPath.toString());
        LOGGER.info("Exon output file: " + exonPath.toString());
        LOGGER.info("Unclassified output file: " + unclassifiedPath.toString());
        LOGGER.info("Padding exons with " + exonPadding + " bases.");
        if (stranded) {
            LOGGER.info("Considering strandedness when calculating overlap.");
        } else {
            LOGGER.info("Not considering strandedness when calculating overlap.");
        }
        
        if (debug) {
            LOGGER.info("Debug mode is on.");
            LOGGER.info("Writing gene-bodies of BED records to " +
                    DEBUG_GENE_BODIES_PATH.toString());
            LOGGER.info("Writing introns of BED records to " +
                    DEBUG_INTRONS_PATH.toString());
        } else {
            LOGGER.info("Debug mode is off.");
        }
        
        return this;
    }
    
    /**
     * Loads records from a BED file into memory.
     * @return this IntronExonSplitter for method chaining
     */
    private IntronExonSplitter loadGenesAndIntrons() {
        loadGenes();
        findIntrons();
        if (debug) writeDebugOutput();
        return this;
    }
    
    /**
     * Loads records from a BED file as both genes and gene-bodies.
     */
    private void loadGenes() {
        LOGGER.info("Loading gene annotations.");
        try (BedParser p = new BedParser(geneBedPath)) {
            p.stream().forEach(x -> {
                if (stranded) {
                    genes.add(x);
                    geneBodies.add(x.getBody());
                } else {
                    genes.add(new Annotation(x, Strand.BOTH));
                    geneBodies.add(new Annotation(x.getBody(), Strand.BOTH));
                }
            });
        }

        genes.stream()
            .forEach(x -> paddedGenes.add(padAnnotation(x, exonPadding)));
        
        LOGGER.info("Finished loading gene annotations.");
    }
    
    /**
     * Calculates the introns from the genes and stores them in memory.
     */
    private void findIntrons() {
        LOGGER.info("Creating intron annotations.");
        paddedGenes.stream()
            .flatMap(x -> x.getIntronStream())
            .forEach(x -> introns.add(x));
        LOGGER.info("Finished creating introns.");
    }

    /**
     * Writes gene bodies and introns to disk as BED files for checking output
     * in IGV.
     */
    private void writeDebugOutput() {
        LOGGER.info("Writing genes and introns for debug check.");
        
        try (BufferedWriter bw = Files.newBufferedWriter(
                DEBUG_INTRONS_PATH, StandardCharsets.US_ASCII)) {
            introns.stream().forEach(x -> {
             try {
                 bw.write(x.toFormattedString());
             } catch (IOException e) {
                 LOGGER.log(Level.WARNING, "Exception caught writing debug output", e);
             }
         });
        } catch (IOException e) {
            LOGGER.log(Level.WARNING, "Exception caught writing debug output", e);
        }
        
        try (BufferedWriter bw = Files.newBufferedWriter(
                DEBUG_GENE_BODIES_PATH, StandardCharsets.US_ASCII)) {
            geneBodies.stream().forEach(x -> {
             try {
                 bw.write(x.toFormattedString());
             } catch (IOException e) {
                 LOGGER.log(Level.WARNING, "Exception caught writing debug output", e);
             }
         });
        } catch (IOException e) {
            LOGGER.log(Level.WARNING, "Exception caught writing debug output", e);
        }
    }
    
    /**
     * Convenience method to get successfully aligned fragments from a BAM file
     * as a stream.
     */
    private Stream<? extends Alignment> getAlignments() {
        LOGGER.info("Parsing BAM file");
        return bamParser.getAlignmentStream();
    }
    
    /**
     * Writes the alignment to the appropriate writer.
     * <p>
     * All logic to assign to the appropriate writer happens here.
     */
    private void writeAlignment(Alignment alignment) {
        
        // Can't assign to a gene. Write as unclassified.
        if (!geneBodies.overlaps(alignment)) {
            unclassifiedWriter.writeSamRecord(alignment);
            return;
        }

        Iterator<Annotated> overlappingGenes = paddedGenes.overlappers(alignment);
        while (overlappingGenes.hasNext()) {
            Annotated overlappingGene = overlappingGenes.next();
            if (overlappingGene.contains(alignment)) {
                if (alignment.isSpliced()) {
                    exonWriter.writeSamRecord(alignment);
                    return;
                }
                
                if (introns.overlaps(alignment)) {
                    unclassifiedWriter.writeSamRecord(alignment);
                    return;
                }
                
                exonWriter.writeSamRecord(alignment);
                return;
            }
        }

        if (introns.overlaps(alignment)) {
            intronWriter.writeSamRecord(alignment);
            return;
        }
            
        unclassifiedWriter.writeSamRecord(alignment);
    }
    
    /**
     * Pads each exon in an annotation on both ends
     * @param annot - the annotation
     * @param padding - the number of bases to pad on each end
     * @return the new padded annotation
     */
    private Annotated padAnnotation(Annotated annot, int padding) {
        AnnotationBuilder builder = Annotation.builder();
        annot.getBlockStream()
            .forEach(x -> builder.addAnnotation(new Annotation(
                    x.getReferenceName(), x.getStart() - exonPadding,
                    x.getEnd() + exonPadding, x.getStrand())));
        return builder.build();
    }
    
    @Override
    public void close() {
        if (intronWriter != null) intronWriter.close();
        if (exonWriter != null) exonWriter.close();
        if (unclassifiedWriter != null) unclassifiedWriter.close();
    }
}