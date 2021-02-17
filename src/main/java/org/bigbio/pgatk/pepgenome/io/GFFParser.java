package org.bigbio.pgatk.pepgenome.io;

import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.PepGenomeTool;
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

    // TODO Edited - Translation start point offset (For use with unannotated peptides)
    private static int translationOffset = 0;


    // GFF3 Patterns
    private static final Pattern GFFIDPATTERN = Pattern.compile("ID=([^;]*)"); // ID tag
    private static final Pattern GFFPARENTPATTERN = Pattern.compile("Parent=([^;]*)"); // Parent tag

    // Singleton get_instance method.
    public static GFFParser get_instance() {
        if (instance == null) {
            instance = new GFFParser();
        }
        return instance;
    }

    // Opens file stream, returns true if successful
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





    // Returns true if in the feature line at position 6 there is a + (plus strand)
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

    // TODO Slight edit ("transcript" --> "mRNA")
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
        //TODO Added remainingProteinLength variable (before parsing loop)
        int remainingProteinLength = 0;

        Coordinates prevProteinCoordinates = new Coordinates();
        Assembly assem = Assembly.none;
        ArrayList<String> tokens;

        // TODO Added hash map of transcript IDs to gene IDs so that gene ID may be retrieved for exons via transcripts (children of genes, parents of exons).
        HashMap<String, String> idMap = new HashMap<>();

        //TODO Parser read loop starts here
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("#")) {
                continue;
            }

            // Convert GFF line into 9 tokens
            tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(line, "\t")));

            // GENE
            if (is_next_gene(tokens)) {
                Assembly assemtemp = mapping.add_gene_from_gtf(line); // check this TODO Note: add_gene_from_gtf used here.
                if (assem == Assembly.none) {
                    if (assemtemp == Assembly.patchhaploscaff) {
                        assem = assemtemp;
                    }
                }
            }

            // TRANSCRIPT
            // String transcriptId = extract_transcript_id(line); // TODO Edited: Removed line.  Inserted replacement a few lines down.****
            String transcriptId = ""; // TODO Edited: Added line to maintain correct scope. ****
            if (is_next_transcript(tokens)) {

                // TODO Edited - Added transcript and associated gene id to idMap for later access.
                String geneId = extract_id(line, GFFPARENTPATTERN);
                transcriptId = extract_id(line, GFFIDPATTERN); //TODO Edited - Replacement for removed line above.****
                idMap.put(transcriptId,geneId); // Places transcript id and gene id into hash map.

                exonID = "";
                mapping.add_transcript_id_to_gene(line); //TODO Note - MAPPEDPEPTIDES add_transcript_id_to_gene used here.
                if (proteinEntry != null) {
                    proteinEntry.set_coordinate_map(coordinatesMap);
                }
                // TODO Count length of of protein, check against exon peptide lengths.
                // TODO Length of protein here
                proteinEntry = coordwrapper.lookup_entry(transcriptId);
                if (proteinEntry == null) {
                    Log.info("ERROR: No entry for transcript ID: " + transcriptId);
                    continue;
                }
                //TODO remainingProteinLength --> Assignment during transcript processing
                remainingProteinLength = proteinEntry.get_sequence().length();

                prevProteinCoordinates = new Coordinates();
                prevProteinCoordinates.setCterm(Offset.off3); //cterm offset (see enum common.Offset)
                prevProteinCoordinates.setNterm(Offset.off3); //nterm offset (see enum common.Offset)
                prevProteinCoordinates.setStart(0); //the start position
                prevProteinCoordinates.setEnd(0); //the end position
                coordinatesMap = new ArrayList<>();

                // TODO Edit - Retrieve transcript's translation offset based on its ID.

                if (PepGenomeTool.m_translation_offset_map.get(transcriptId) == null) {
                    translationOffset = 0;
                }
                else {
                    translationOffset = PepGenomeTool.m_translation_offset_map.get(transcriptId);
                }

                // For testing:
                //translationOffset = PepGenomeTool.m_offsetMap.get("alt_3prime.4188_iso2"); // 0 (Passed)
                //translationOffset = PepGenomeTool.m_translation_offset_map.get("alt_3prime.4201_iso2"); // 385 (Passed)


                // Offset now set before exons are processed.  Processing am exon involves examining its length in nucleotides.
                // This is compared to offset value to determine whether exon is deleted, changed or left alone.


                // EXON
            } else if (is_exon(tokens)) {

                if (PepGenomeTool.useExonCoords) {

                    // Altered CDS block

                    // TODO Finish implementing the following (Finished, now check everything is correct, and test):
                    // Chris: Change genCoord based on JAMES's offset from translation (e.g.
                    // exon completely untranslated (James offset longer than exon) -> remove -> nothing else required to do.
                    // exon partially translated (James offset shorter than exon) -> subtract offset from genCoord.start / add to genCoord.end (depending on strand).
                    // change offset -> James offset - length of exon or length of untranslated part of exon (ie until James offset is 0).
                    // This means change the offset as parts are removed and subtracted from the  previous exons.  Therefore J's offset will need to be updated upon analysing a given exon.

                    /*
                    Directions from Chris / Explanations for following code edits:
                    Translation offset from FASTA headers to be applied here.

                    Transcript read by GFF Parser, transcript ID searched for within translation offset map (PepGenomeTool.java).
                    If map entry exists for transcript ID, translation offset read.

                    Assuming annotation file is arranged correctly, exons will follow.  Check if this will always be the case.

                    So while tracking the offset value for the given transcript, checking each exon has a parent matching the trans ID...

                    If exon totally untranslated:
                    (ie the offset is greater than the length of the exon in )
                        Remove exon
                    Else if exon partially untranslated:
                        Exon remains but with genomic coordinates adjusted by subtracting translation offset.
                        Adjust offset value accordingly before moving onto next exon.


                    Also mind strand direction (Forward, Reverse).
                     */

                    // TODO Call 2 Notes
                    // Generate input files for every case (10):
                    // ...
                    // ...Additional exon that isnt translated
                    // Same for reverse
                    // Extreme cases


                    // TODO Edited - Check exon has parent ID matching  the last transcript ID - Ensure offset is being applied correctly.  Possibly unnecessary.
                    if (extract_id(line, GFFPARENTPATTERN).equals(transcriptId)) {

                        // Extract genomic coordinates from exon line
                        GenomeCoordinates genomeCoordinates = Utils.extract_coordinates_from_gtf_line(tokens);

                        // TODO Edited - Determine length between start and end points of exon (Moved position of this code, original position commented out below.)
                        //  Determining length of exon
                        int length = 0;
                        if (is_first_strand(tokens)) {
                            length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1;
                        } else if (!is_first_strand(tokens)) {
                            length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1;
                        }

                        // TODO Edit - Check length of exon against offset.  Adjust coords accordingly.
                        //TODO Edit - Check remaining protein length, 0 or below indicates all following exons are untranslated.

                        // UNTRANSLATED EXON
                        if (translationOffset > length || remainingProteinLength <= 0 ) {
                            // Offset value greater than total exon length - Exon  is completely untranslated

                            // Adjust offset
                            translationOffset = translationOffset - length;
                            System.out.println(translationOffset);

                            continue;


                            // PARTIALLY TRANSLATED EXON
                        } else if (translationOffset > 0) {
                            // Offset value not greater than exon length, but still greater than 0 - Exon is partially translated.
                            // TODO Adjust genomic coords - Checked
                            if (is_first_strand(tokens)) {
                                genomeCoordinates.setStart(genomeCoordinates.getStart() + translationOffset);
                            } else {
                                genomeCoordinates.setEnd(genomeCoordinates.getEnd() - translationOffset);
                            }

                            // Adjust offset
                            translationOffset = 0;
                        }

                        // PARTIALLY/FULLY TRANSLATED EXON
                        genomeCoordinates.setTranscriptid(transcriptId);  // Using previous transcript's ID
                        exonID = extract_exon_id(line);  // Extracted from current exon line
                        genomeCoordinates.setExonid(exonID);
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

                        //  Determining length (Original position, moved above so length can be used with translation offset.
                    /*
                    int length = 0;
                    if (is_first_strand(tokens)) {
                        length = genomeCoordinates.getEnd() - genomeCoordinates.getStart() + 1;
                    } else if (!is_first_strand(tokens)) {
                        length = genomeCoordinates.getStart() - genomeCoordinates.getEnd() + 1;
                    }

                     */

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

                        // TODO Change remaining protein length (RPL) to subtract the length of the peptide being processed.
                        remainingProteinLength = remainingProteinLength - peplength;



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

                else {
                    // System.out.println("Annotation: Using CDS Coords");
                    // CDS analysed separately in own branch.
                    // TODO Process exon as in original parser.
                }



        } else if (is_cds(tokens)) {

                // Normal CDS Block, edited to work with GFF3 but no further changes.
            GenomeCoordinates genomeCoordinates = Utils.extract_coordinates_from_gtf_line(tokens); // Chris - Should be fine
            genomeCoordinates.setTranscriptid(transcriptId);
            String tmp_exonID = extract_id(line,GFFPARENTPATTERN); // TODO Edited to use Parent ID rather than own ID (Previously extract_exon_id which wont work in the CDS branch.)

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

    // TODO Edited -looks for the text specified and returns the ID. (Working)
    public static String extract_gene_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GFFIDPATTERN);
    }

    // TODO Edited -looks for the text specified and returns the ID. (Working)
    public static String extract_transcript_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GFFIDPATTERN);
    }

    //TODO Edited -looks for the text specified and returns the ID. (Working)
    public static String extract_exon_id(String gtfGeneLine) {
        return extract_id(gtfGeneLine, GFFIDPATTERN);
    }

    // TODO Edited - General extract_id method.  Other extract ID methods will now call this. (Working)
    public static String extract_id(String gtfGeneLine, Pattern pattern) {
        String value = "";
        Matcher matcher = pattern.matcher(gtfGeneLine);
        if (matcher.find()) {
            value = matcher.group(1);
        }
        return value;
    }

    // TODO Edited - extract_gene_name method moved from GeneEntry to GFF Parser and edited.  Original version moved to GTF Parser.
    //extracts the gene symbol
    public static String extract_gene_name(List<String> tokens) {
        String value = "";
        if (tokens.size() >= 9) {
            List<String> res = GeneEntry.extract_by_tag("gene_name", tokens.get(8));

            if (res.size() == 1) {
                value = res.get(0);
            }
        }
        return value;
    }

}
