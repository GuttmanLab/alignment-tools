package edu.caltech.lncrna.alignment.programs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.Optional;
import java.util.logging.Logger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.caltech.lncrna.bio.annotation.Annotated;
import edu.caltech.lncrna.bio.annotation.Annotation;
import edu.caltech.lncrna.bio.annotation.Annotation.AnnotationBuilder;
import edu.caltech.lncrna.bio.annotation.BedFileRecord;
import edu.caltech.lncrna.bio.annotation.Populated;
import edu.caltech.lncrna.bio.annotation.Strand;
import edu.caltech.lncrna.bio.annotation.WindowIterator;
import edu.caltech.lncrna.bio.io.BedParser;

/**
 * This program accepts a BED-formatted mask file as input and outputs
 * (for default parameters) a BED6 file of 1-Mb tiles for which the score is
 * equal to the percentage of the tile that was masked.
 * <p>
 * For example, a line might be:
 * <p>
 * <code>chr1  5000000  6000000  .  0.9500  .</code>
 * <p>
 * if the chr1:5000000-6000000 tile is 95% masked.
 */
public final class MaskQuantifier {

    private final Path maskPath;
    private final Path outputPath;
    private final int windowSize;
    private final int staggerSize;
    
    private static final Logger LOGGER = Logger.getLogger("MaskQuantifier");
    public static final String VERSION = "1.0.0";
    private static final int DEFAULT_WINDOW_SIZE = 1000000; // 1 Mb
    private static final int NUMBER_BED_OUTPUT_FIELDS = 6;
    
    public static void main(String[] args) throws IOException {
        long startTime = System.currentTimeMillis();
        CommandLine cmd = parseArgs(args);

        new MaskQuantifier(cmd).countOverWindows();
        
        LOGGER.info("Program complete.");
        LOGGER.info((System.currentTimeMillis() - startTime) + " milliseconds elapsed.");
    }
    
    public MaskQuantifier(CommandLine cmd) {
        maskPath = Paths.get(cmd.getOptionValue("mask"));
        outputPath = Paths.get(cmd.getOptionValue("output"));
        
        windowSize = cmd.hasOption("window")
                ? Integer.parseInt(cmd.getOptionValue("window"))
                : DEFAULT_WINDOW_SIZE;
                
        staggerSize = cmd.hasOption("stagger")
                ? Integer.parseInt(cmd.getOptionValue("stagger"))
                : windowSize;
    }
    
    private static CommandLine parseArgs(String[] args) {
        
        Option versionOption = Option.builder("v")
                .longOpt("version")
                .desc("show version information")
                .hasArg(false)
                .required(false)
                .build();
        
        Option helpOption = Option.builder("h")
                .longOpt("help")
                .desc("show usage information")
                .hasArg(false)
                .required(false)
                .build();
        
        Option maskOption = Option.builder()
                .longOpt("mask")
                .desc("mask BED file")
                .hasArg()
                .required()
                .build();

        Option outputOption = Option.builder()
                .longOpt("output")
                .desc("output BED file")
                .hasArg()
                .required()
                .build();
        
        Option windowOption = Option.builder()
                .longOpt("window")
                .desc("length of windows to tile genome (defaults to 1 Mb)")
                .argName("NUM")
                .hasArg()
                .required(false)
                .build();
        
        Option staggerOption = Option.builder()
                .longOpt("stagger")
                .desc("tiling offset (defaults to window-length for 1x " +
                        "coverage)")
                .argName("NUM")
                .hasArg()
                .required(false)
                .build();
        
        Options helpOptions = new Options().addOption(helpOption);
        Options versionOptions = new Options().addOption(versionOption);
        
        Options mainOptions = new Options()
                .addOption(maskOption)
                .addOption(outputOption)
                .addOption(windowOption)
                .addOption(staggerOption);
        
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
                formatter.printHelp("test", allOptions);
                System.exit(0);
            }
            cmds = parser.parse(versionOptions, args, true);
            if (cmds.getOptions().length == 1) {
                System.out.println(VERSION);
                System.exit(0);
            }
            rtrn = parser.parse(mainOptions, args);
        } catch (ParseException e) {
            formatter.printHelp("test", allOptions);
            System.exit(1);
        }
        
        if (rtrn == null) {
            LOGGER.severe("An unknown error occurred while parsing command line arguments");
            System.exit(1);
        }
        
        return rtrn;
    }
    
    private void countOverWindows() throws IOException {
        try (BedParser bp = new BedParser(maskPath);
             BufferedWriter bw = Files.newBufferedWriter(outputPath)) {
            
            WindowIterator<BedFileRecord> windows =
                    new WindowIterator<>(bp, windowSize, staggerSize);
            
            while (windows.hasNext()) {
                Populated<BedFileRecord> window = windows.next();

                AnnotationBuilder builder = new AnnotationBuilder();
                Iterator<BedFileRecord> pop = window.getPopulation();
                while (pop.hasNext()) {
                    builder.addAnnotation(new Annotation(pop.next(), Strand.BOTH));
                }
                Annotated mask = builder.build();
                
                Optional<Annotated> unmasked = window.minus(mask);
                float numberUnmaskedBases = unmasked.isPresent()
                        ? unmasked.get().getSize()
                        : 0;
                
                float percentMasked = 1 - numberUnmaskedBases / window.getSpan();
                String name = window.getReferenceName() + ":" + window.getStart() +
                        "-" + window.getEnd();
                bw.write(BedFileRecord.builder()
                        .addAnnotation(window)
                        .addScore(percentMasked)
                        .addName(name)
                        .build()
                        .toFormattedBedString(NUMBER_BED_OUTPUT_FIELDS));
            }
        }
    }
}
