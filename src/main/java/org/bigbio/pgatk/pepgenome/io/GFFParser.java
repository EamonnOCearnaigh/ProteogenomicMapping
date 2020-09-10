package org.bigbio.pgatk.pepgenome.io;

import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.common.*;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.mortbay.log.Log;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// GFF3 Parser
public class GFFParser {


    private static Logger log = LoggerFactory.getLogger(GTFParser.class);

    // Input stream
    private BufferedReader reader;

    private FileInputStream ifs;

    // Current line
    private String line;

    public static GFFParser instance;

    private GFFParser() {
    }


    // GFF3 - All relevant feature types are reduced to the tags 'ID' and 'Parent' rather than GTF's 'gene ID' and 'transcript ID'.
    private static Pattern GFFIDPATTERN = Pattern.compile("ID=([^;.]*)"); // ID tag
    private static Pattern GFFPARENTPATTERN = Pattern.compile("Parent=([^;.]*)"); // Parent tag

    // Singleton get_instance method.
    public static GFFParser get_instance() {
        if (instance == null) {
            instance = new GFFParser();
        }
        return instance;
    }

    // Opens file stream, returns true if successful - Same
    private boolean open(String file) throws Exception {
        if (reader == null) {
            line = "";
            ifs = new FileInputStream(file);
            reader = new BufferedReader(new InputStreamReader(ifs));
        }
        boolean status = true;
        try{
            status = reader.ready();
        }catch (IOException ex){
            log.debug("The gff3 file stream is closed -- " + file);
            ifs = new FileInputStream(file);
            reader = new BufferedReader(new InputStreamReader(ifs));
        }
        return status;
    }

    // Closes the file stream
    private void close() throws Exception {
        line = "";
        reader.close();
        ifs.close();
    }





    // Returns true if in the GTF at position 6 there is a + (plus strand) - Same as GTF
    private static boolean is_first_strand(List<String> tokens) {
        return tokens.get(6).equals("+");
    }

    // Returns true if line feature type is CDS
    private static boolean is_cds(List<String> tokens) {
        return tokens.get(2).equals("CDS");
    }

    // Returns true if line feature type is exon
    private static boolean is_exon(List<String> tokens) {
        return tokens.get(2).equals("exon");
    }

    // Returns true if line feature type is transcript - Switched "transcript" to "mRNA"
    private static boolean is_next_transcript(List<String> tokens) { return tokens.get(2).equals("mRNA"); }

    // Returns true if line feature type is gene
    private static boolean is_next_gene(List<String> tokens) {
        return tokens.get(2).equals("gene");
    }


    // Reads a gff3 file and parses it into CoordinateWrapper and MappedPeptides.

    public final Assembly read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping) throws Exception {
        if (!open(file)) {
            throw new IllegalStateException("Problem in reading GFF3 file");
        }

        String exonID = "";
        ProteinEntry proteinEntry = null;
        ArrayList<Tuple<Coordinates, GenomeCoordinates>> coordinatesMap = new ArrayList<>();

        Coordinates prevProteinCoordinates = new Coordinates();
        Assembly assem = Assembly.none;
        ArrayList<String> tokens;

        // TODO Added hash map of transcript IDs to gene IDs so that gene ID may be retrieved for exons across two generations.
        HashMap<String, String> idMap = new HashMap<String, String>();

        while ((line = reader.readLine()) != null) {
            if (line.startsWith("#")) {
                continue;
            }

            // Convert GFF line into 9 tokens
            tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(line, "\t")));

            // GENE
            if (is_next_gene(tokens)) {
                Assembly assemtemp = mapping.add_gene_from_gtf(line); // check this TODO Note - add_gene_from_gtf used here.
                if (assem == Assembly.none) {
                    if (assemtemp == Assembly.patchhaploscaff) {
                        assem = assemtemp;
                    }
                }
            }

            // TRANSCRIPT
            String transcriptId = extract_transcript_id(line); //TODO Note - GENEENTRY extract_transcript_id_gff used here.
            if (is_next_transcript(tokens)) {

                // TODO Edit - Add transcript and associated gene to idMap for later access.
                String geneId = extract_id(line, GFFPARENTPATTERN);
                idMap.put(transcriptId,geneId);

                exonID = "";
                mapping.add_transcript_id_to_gene(line); //TODO Note - MAPPEDPEPTIDES add_transcript_id_to_gene used here.
                if (proteinEntry != null) {
                    proteinEntry.set_coordinate_map(coordinatesMap);
                }


                proteinEntry = coordwrapper.lookup_entry(transcriptId);
                if (proteinEntry == null) {
                    Log.info("ERROR: No entry for transcript ID: " + transcriptId);
                    continue;
                }

                prevProteinCoordinates = new Coordinates();
                prevProteinCoordinates.setCterm(Offset.off3); //cterm offset (see enum common.Offset)
                prevProteinCoordinates.setNterm(Offset.off3); //nterm offset (see enum common.Offset)
                prevProteinCoordinates.setStart(0); //the start position
                prevProteinCoordinates.setEnd(0); //the end position
                coordinatesMap = new ArrayList<>();


                // EXON
            } else if (is_exon(tokens)) {

                //TODO Removed line: exonID = extract_exon_id(line);  Replaced with former CDS block:


                // Sometimes there will not be CDS, so we will need to do all of the work in Exon branch and use exon line to get gene coords.  But only for first couple of exons within a transcript y
                // special part using the output of james, looking at start of transcript, where does translation start.
                // Manually change gene coords.

                // Moved CDS block here.
                GenomeCoordinates genomeCoordinates = Utils.extract_coordinates_from_gtf_line(tokens); // Chris - Should be fine
                // TODO Next point of focus *****************
                // Chris: Change genCoord based on JAMES's offset from translation (e.g.
                // exon completely untranslated (James offset longer than exon) -> remove -> nothing else required to do.
                // exon partially translated (James offset shorter than exon) -> subtract offset from genCoord.start / add to genCoord.end (depending on strand).
                // change offset -> James offset - length of exon or length of untranslated part of exon (ie until James offset is 0).
                // This means change the offset as parts are removed and subtracted from the  previous exons.  Therefore J's offset will need to be updated.

                genomeCoordinates.setTranscriptid(transcriptId);
                String tmp_exonID = extract_exon_id(line);

                if (tmp_exonID.equals("")) {
                    tmp_exonID = exonID;
                }

                genomeCoordinates.setExonid(tmp_exonID);
                Coordinates proteinCoordinates = new Coordinates();

                // Get N term from previous exon
                if (genomeCoordinates.getFrame() != Frame.unknown) {
                    proteinCoordinates.setNterm(Offset.forValue(genomeCoordinates.getFrame().getValue()));
                } else {
                    if (prevProteinCoordinates.getCterm() != Offset.off3) {
                        proteinCoordinates.setNterm(Offset.forValue(3 - prevProteinCoordinates.getCterm().getValue()));
                    } else {
                        proteinCoordinates.setNterm(Offset.off3);
                    }
                }

                //  Determining length
                int length = 0;
                if (is_first_strand(tokens)) {
                    length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1;
                } else if (!is_first_strand(tokens)) {
                    length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1;
                }

                // Calculating C term
                if (length % 3 == 0) {
                    if (proteinCoordinates.getNterm() != Offset.off3) {
                        proteinCoordinates.setCterm(Offset.forValue(3 - proteinCoordinates.getNterm().getValue()));
                    } else {
                        proteinCoordinates.setCterm(Offset.off3);
                    }
                } else if (length % 3 == 2) {
                    if (proteinCoordinates.getNterm() == Offset.off3) {
                        proteinCoordinates.setCterm(Offset.off2);
                    } else if (proteinCoordinates.getNterm() == Offset.off2) {
                        proteinCoordinates.setCterm(Offset.off3);
                    } else if (proteinCoordinates.getNterm() == Offset.off1) {
                        proteinCoordinates.setCterm(Offset.off1);
                    }
                } else if (length % 3 == 1) {
                    if (proteinCoordinates.getNterm() == Offset.off3) {
                        proteinCoordinates.setCterm(Offset.off1);
                    } else if (proteinCoordinates.getNterm() == Offset.off1) {
                        proteinCoordinates.setCterm(Offset.off3);
                    } else if (proteinCoordinates.getNterm() == Offset.off2) {
                        proteinCoordinates.setCterm(Offset.off2);
                    }
                }

                // Calculate protein coordinates
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    proteinCoordinates.setStart(prevProteinCoordinates.getEnd());
                } else {
                    if (prevProteinCoordinates.getEnd() == 0 && coordinatesMap.isEmpty()) {
                        proteinCoordinates.setStart(0);
                    } else {
                        proteinCoordinates.setStart(prevProteinCoordinates.getEnd() + 1);
                    }
                }

                int offsets = 0;
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    offsets = offsets + proteinCoordinates.getNterm().getValue();
                }

                if (is_first_strand(tokens)) {
                    length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1 - offsets;
                } else if (!is_first_strand(tokens)) {
                    length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1 - offsets;
                }

                int peplength = length / 3;

                int pepend = proteinCoordinates.getStart() + peplength - 1;
                if (proteinCoordinates.getCterm() != Offset.off3) {
                    ++pepend;
                }
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    ++pepend;
                }

                proteinCoordinates.setEnd(pepend);

                prevProteinCoordinates = proteinCoordinates;

                coordinatesMap.add(new Tuple<>(proteinCoordinates, genomeCoordinates));

            }

        }

        if (proteinEntry != null) {
            proteinEntry.set_coordinate_map(coordinatesMap);
        }
        close();
        return assem;
    }

        /*

             // CDS
            } else if (is_cds(tokens)) {
                GenomeCoordinates genomeCoordinates = Utils.extract_coordinates_from_gtf_line(tokens); // Chris - Should be fine
                genomeCoordinates.setTranscriptid(transcriptId);
                String tmp_exonID = extract_exon_id(line);

                if (tmp_exonID.equals("")){
                    tmp_exonID = exonID;
                }

                genomeCoordinates.setExonid(tmp_exonID);
                Coordinates proteinCoordinates = new Coordinates();

                // Get N term from previous exon
                if (genomeCoordinates.getFrame() != Frame.unknown) {
                    proteinCoordinates.setNterm(Offset.forValue(genomeCoordinates.getFrame().getValue()));
                } else {
                    if (prevProteinCoordinates.getCterm() != Offset.off3) {
                        proteinCoordinates.setNterm(Offset.forValue(3 - prevProteinCoordinates.getCterm().getValue()));
                    } else {
                        proteinCoordinates.setNterm(Offset.off3);
                    }
                }

                //  Determining length
                int length = 0;
                if (is_first_strand(tokens)) {
                    length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1;
                } else if (!is_first_strand(tokens)) {
                    length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1;
                }

                // Calculating C term
                if (length % 3 == 0) {
                    if (proteinCoordinates.getNterm() != Offset.off3) {
                        proteinCoordinates.setCterm(Offset.forValue(3 - proteinCoordinates.getNterm().getValue()));
                    } else {
                        proteinCoordinates.setCterm(Offset.off3);
                    }
                } else if (length % 3 == 2) {
                    if (proteinCoordinates.getNterm() == Offset.off3) {
                        proteinCoordinates.setCterm(Offset.off2);
                    } else if (proteinCoordinates.getNterm() == Offset.off2) {
                        proteinCoordinates.setCterm(Offset.off3);
                    } else if (proteinCoordinates.getNterm() == Offset.off1) {
                        proteinCoordinates.setCterm(Offset.off1);
                    }
                } else if (length % 3 == 1) {
                    if (proteinCoordinates.getNterm() == Offset.off3) {
                        proteinCoordinates.setCterm(Offset.off1);
                    } else if (proteinCoordinates.getNterm() == Offset.off1) {
                        proteinCoordinates.setCterm(Offset.off3);
                    } else if (proteinCoordinates.getNterm() == Offset.off2) {
                        proteinCoordinates.setCterm(Offset.off2);
                    }
                }

                // Calculate protein coordinates
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    proteinCoordinates.setStart(prevProteinCoordinates.getEnd());
                } else {
                    if (prevProteinCoordinates.getEnd() == 0 && coordinatesMap.isEmpty()) {
                        proteinCoordinates.setStart(0);
                    } else {
                        proteinCoordinates.setStart(prevProteinCoordinates.getEnd() + 1);
                    }
                }

                int offsets = 0;
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    offsets = offsets + proteinCoordinates.getNterm().getValue();
                }

                if (is_first_strand(tokens)) {
                    length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1 - offsets;
                } else if (!is_first_strand(tokens)) {
                    length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1 - offsets;
                }

                int peplength = length / 3;

                int pepend = proteinCoordinates.getStart() + peplength - 1;
                if (proteinCoordinates.getCterm() != Offset.off3) {
                    ++pepend;
                }
                if (proteinCoordinates.getNterm() != Offset.off3) {
                    ++pepend;
                }

                proteinCoordinates.setEnd(pepend);

                prevProteinCoordinates = proteinCoordinates;

                coordinatesMap.add(new Tuple<>(proteinCoordinates, genomeCoordinates));

            }

        }

        if (proteinEntry != null) {
            proteinEntry.set_coordinate_map(coordinatesMap);
        }
        close();
        return assem;
    }

         */


    //TODO Edited -looks for the text specified GENEPATTERN and returns the ID.
    public static String extract_gene_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GFFIDPATTERN);
    }

    //TODO Edited -looks for the text specified in TRNASCRIPTPATTERN and returns the ID.
    public static String extract_transcript_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GFFIDPATTERN);
    }

    //TODO Edited -looks for the text specified in EXONPATTERN and returns the ID.
    public static String extract_exon_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GFFIDPATTERN);
    }

    // TODO Edited - General extract_id method.  Other extract ID methods will now call this.
    public static String extract_id(String gtfGeneLine, Pattern pattern) {
        String value = "";
        Matcher matcher = pattern.matcher(gtfGeneLine);
        if (matcher.find()) {
            value = matcher.group(1);
        }
        return value;
    }

}
