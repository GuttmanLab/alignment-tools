package edu.caltech.lncrna.alignment.programs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Logger;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.caltech.lncrna.bio.alignment.Aligned;
import edu.caltech.lncrna.bio.alignment.Alignment;
import edu.caltech.lncrna.bio.annotation.WindowIterator;
import edu.caltech.lncrna.bio.io.BamParser;

public class WindowCounter {

    private final Path inputPath;
    private final Path outputPath;
    private final int windowSize;
    private final int staggerSize;
    
    private static final Logger LOGGER = Logger.getLogger("WindowCounter");
    public static final String VERSION = "1.0.0";
    private static final int DEFAULT_WINDOW_SIZE = 1000000; // 1 Mb
    private static final int NUMBER_BED_OUTPUT_FIELDS = 6;
    
    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        CommandLine cmd = parseArgs(args);
        new WindowCounter(cmd).getWindowsAndPrint();
        LOGGER.info("Program complete.");
        LOGGER.info((System.currentTimeMillis() - startTime) + " milliseconds elapsed.");
    }
    
    public WindowCounter(CommandLine cmd) {
        inputPath = Paths.get(cmd.getOptionValue("input"));
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
                .longOpt("input")
                .desc("input BAM file")
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

    private void getWindowsAndPrint() {
        try (BamParser<? extends Aligned<? extends Alignment>> bp = 
                BamParser.newInstance(inputPath);
             BufferedWriter bw = Files.newBufferedWriter(outputPath)) {
            
            WindowIterator<? extends Alignment> windows =
                    new WindowIterator<>(bp.getAlignmentIterator(), windowSize,
                    staggerSize);
            
            while (windows.hasNext()) {
                bw.write(windows.next().toFormattedBedString(NUMBER_BED_OUTPUT_FIELDS));
            }
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
}
