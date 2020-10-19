package org.bigbio.pgatk.pepgenome.common;

import org.bigbio.pgatk.pepgenome.PepGenomeTool;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
// TODO Edited- Incorrect patterns and gene/transcript extraction methods  - were not recognising James's FASTA headers as correct input.  Commented out originals and altered copies, etc.
public class ProteinEntry implements Serializable {
    private static final long serialVersionUID = -1732196455282495216L;
    //whole fasta header
    private String m_fasta_header;
    //transcript id
    private String m_transcript_id;
    //gene id
    private String m_gene_id;
    // TODO Edited - Added offset field
    private int m_offset = 0;
    //the AA sequence.
    private String m_aa_sequence;
    // Original regex patterns
    private static Pattern GENEPATTERN = Pattern.compile("gene:([^\\s\\.]*)[^\\s]*\\s"); // Gene ID tag
    private static Pattern TRANSCRIPTPATTERN = Pattern.compile("transcript:([^\\s\\.]*)[^\\s]*\\s"); // Transcript ID tag

    // TODO Edited - Spladder regex patterns
    private static Pattern spladderGENEPATTERN = Pattern.compile("gene=(.*)");
    private static Pattern spladderTRANSCRIPTTPATTERN = Pattern.compile(">(.*)\\|");
    private static Pattern spladderOFFSETPATTERN = Pattern.compile(" offset=(\\d*) ");

    //std::multimap <Coordinates (protein coordinates), GenomeCoordinates(corresponding genomic coordinates), Coordinates (passing this as third argument will use the Coordinates::operator() as comparator)>
    //the first Coordinate are the coordinates of exons within the protein and the GenomeCoordinate is its corresponding location in the genome.
    private ArrayList<Tuple<Coordinates, GenomeCoordinates>> m_coordinates_map;
    //check if the coding sequence is dividable by 3bp and not offset due to incomplete transcript annotation.
    private int m_cds_annotation_correct;

    public ProteinEntry() {
        this.m_fasta_header = "";
        this.m_transcript_id = "";
        this.m_gene_id = "";
        this.m_aa_sequence = "";
        this.m_coordinates_map = new ArrayList<>();
        this.m_cds_annotation_correct = 0;
        this.m_offset = 0;
    }

    public ProteinEntry(String fastaHeader, String AAsequence) {
        init(fastaHeader, AAsequence);
    }

    public ProteinEntry(FastaEntry fastaEntry) {
        init(fastaEntry);
    }

    //QOL: delegates to void init(std::string fastaHeader, std::string AAsequence);
    private void init(FastaEntry fastaEntry) {
        init(fastaEntry.get_header(), fastaEntry.get_sequence());
    }

    //sets all crucial values.
    private void init(String fastaHeader, String AAsequence) {
        if (fastaHeader.substring(0, 1).equals(">")) {
            m_fasta_header = fastaHeader;
            m_transcript_id = extract_transcript_id_fasta(fastaHeader);
            m_gene_id = extract_gene_id_fasta(fastaHeader);
            m_aa_sequence = AAsequence;
            m_coordinates_map = new ArrayList<>();
            m_cds_annotation_correct = 0;
            m_offset = extract_offset_fasta(fastaHeader);
        }
    }

    // TODO Edited method for spladder format compatibility.  Original method preserved below.

    // TODO Edited -  adjusted this method as transcripts were not being read.  Change pattern. (FIXED)
    // gets the transcriptId from a fasta header
    private String extract_transcript_id_fasta(String str) {
    	String value = "";
    	Matcher transcriptMatcher = spladderTRANSCRIPTTPATTERN.matcher(str);

    	// TODO Note: Problem extracting trans ID (FIXED)***

    	if (transcriptMatcher.find()) {
    		value = transcriptMatcher.group(1);
            //System.out.println("Transcript extracted correctly");
            //System.out.println(value);
    	} else {
            System.out.println("Transcript not extracted correctly");
    	    // Temporary measure until pattern fixed.
    	    //value = value.substring(1,21);
    	}
        return value;
    }


    // TODO Edited method for spladder format compatibility.  Original method preserved below.
    //gets the gene id from a fasta header
    private String extract_gene_id_fasta(String str) {
    	String value = "";
        // TODO Original: Matcher genetMatcher = GENEPATTERN.matcher(str);
    	Matcher geneMatcher = spladderGENEPATTERN.matcher(str);
    	if (geneMatcher.find()) {
    		value = geneMatcher.group(1);
    	} else {
            String[] split = str.split(" ");
            value = split[2];
    	}
        return value;
    }

    // TODO Edited - Offset extraction (Working)
    //gets the coding start point offset from a fasta header
    private Integer extract_offset_fasta(String str) {
        String value = "";
        Matcher offsetMatcher = spladderOFFSETPATTERN.matcher(str);
        if (offsetMatcher.find()) {
            //System.out.println("Offset extracted correctly");
            value = offsetMatcher.group(1);
        }

        //TODO Edited - Put transcript ID and offset value into map in PepGenomeTool
        m_offset = Integer.parseInt(value);
        //System.out.println("PROTEIN ENTRY: Adding transcript ID and offset to m_offsetMap");
        //System.out.println(m_transcript_id);
        //System.out.println(m_offset);
        PepGenomeTool.m_offsetMap.put(m_transcript_id, m_offset);

        return m_offset;
    }


    // TODO Original extract_transcript_id_fasta method
    /*
    // gets the transcriptId from a fasta header
    private String extract_transcript_id_fasta(String str) {
        String value = "";
        Matcher transcriptMatcher = TRANSCRIPTPATTERN.matcher(str);
        if (transcriptMatcher.find()) {
            value = transcriptMatcher.group(1);
        } else {
            String[] split = str.split("\\|");
            if(split.length==8) {
            	String[] dotsplit = split[1].split("\\.");
            value = dotsplit[0];
            }
        }
        return value;
    }

     */


    // TODO Original extract_gene_id_fasta method
    /*
    // gets the gene id from a fasta header
    private String extract_gene_id_fasta(String str) {
        String value = "";
        Matcher geneMatcher = GENEPATTERN.matcher(str);
        if (geneMatcher.find()) {
            value = geneMatcher.group(1);
        } else {
            String[] split = str.split("\\|");
            if(split.length==8) {
            String[] dotsplit = split[2].split("\\.");
            value = dotsplit[0];
            }
        }
        return value;
    }

     */


    //returns the transcript_id number of the current protein
    public String get_transcript_id() {
        return m_transcript_id;
    }

    //returns the gene_id number of the current protein
    public String get_gene_id() {
        return m_gene_id;
    }

    //returns the sequence. the sequences are iso-sequences (I and L are converted to J)
    public String get_sequence() {
        return m_aa_sequence;
    }

    //setter for the coordinatesMap
    public void set_coordinate_map(ArrayList<Tuple<Coordinates, GenomeCoordinates>> coordinatesMap) {
        m_coordinates_map = coordinatesMap;
    }

    //getter for CDS_annotation_correct.
    public int get_cds_annotation_correct() {
        return m_cds_annotation_correct;
    }

    //returns the genomic coordinates.
    //mapping function. takes the positions calculated in the KmereMap and generates the genomic coordinates for all peptides of this protein.
    public ArrayList<ArrayList<GenomeCoordinates>> find_coordinates(int peptideseqSize, ArrayList<PositionMismatchT> positions) {
        ArrayList<ArrayList<GenomeCoordinates>> foundCoordinates = new ArrayList<>();
        Coordinates peptideCoordinates = new Coordinates();
        peptideCoordinates.setCterm(Offset.off3);
        peptideCoordinates.setNterm(Offset.off3);

        //iterate all found positions
        for (PositionMismatchT current : positions) {
            peptideCoordinates.setStart(current.position_in_protein());
            peptideCoordinates.setEnd(current.position_in_protein() + (peptideseqSize - 1));
            ArrayList<GenomeCoordinates> single = new ArrayList<>();
            m_coordinates_map.stream().filter(e -> e.getKey().equals(peptideCoordinates))
                    .forEach(fe -> {
                        Tuple<Coordinates, GenomeCoordinates> coordinatesPartial = Utils.get_coordinates(fe.getKey(), fe.getValue(), peptideCoordinates);
                        single.add(coordinatesPartial.getValue());
                    });
            //and has to be done several times to find all peptides.
            foundCoordinates.add(single);
        }
        return foundCoordinates;
    }

    @Override
    public String toString() {
        return "ProteinEntry{" +
                "m_fasta_header='" + m_fasta_header + '\'' +
                ", m_transcript_id='" + m_transcript_id + '\'' +
                ", m_gene_id='" + m_gene_id + '\'' +
                ", m_aa_sequence='" + m_aa_sequence + '\'' +
                ", m_coordinates_map=" + m_coordinates_map +
                ", m_cds_annotation_correct=" + m_cds_annotation_correct +
                '}';
    }
}