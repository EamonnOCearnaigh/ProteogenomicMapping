package org.bigbio.pgatk.pepgenome.io;

import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SparkSession;
import org.bigbio.pgatk.pepgenome.CoordinateWrapper;
import org.bigbio.pgatk.pepgenome.PepGenomeTool;
import org.bigbio.pgatk.pepgenome.common.PeptideEntry;
import org.bigbio.pgatk.pepgenome.common.SparkConfig;
import org.bigbio.pgatk.pepgenome.common.TranscriptsT;
import org.bigbio.pgatk.pepgenome.common.Utils;
import org.bigbio.pgatk.pepgenome.common.maps.MappedPeptides;
import org.bigbio.pgatk.pepgenome.kmer.IKmerMap;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * This class parses peptide input files from Tab delimited files. This classes only read the tab delimited files define by
 * Sample | Peptide | PSMs | Quant
 *
 * @author ypriverol
 */

public class TabInputPeptideFileParser implements PeptideInputReader, Serializable {

    // TODO ||Read method||
    public void read(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {
        if (SparkConfig.getInstance().getMaster() != null) {
            System.out.println("TabInputPeptideFileParser.read(...) - sparkRead(...) method used");
            sparkRead(file, coordwrapper, mapping, unmappedoutput, k);

        } else {
            System.out.println("TabInputPeptideFileParser.read(...) - normalRead(...) method used");
            normalRead(file, coordwrapper, mapping, unmappedoutput, k);

        }

    }

    // TODO ||normalRead method||
    //read function. this reads the peptides input and sets the wheels in motion.
    //this function will set the wheels in motion to find the peptides in the proteins.
    private void normalRead(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {

        // TODO ||Peptide file input||
        FileInputStream ifs = new FileInputStream(file);
        BufferedReader reader = new BufferedReader(new InputStreamReader(ifs));
        FileOutputStream ofs = new FileOutputStream(unmappedoutput);

        String peptide_string;
        String tissue;
        String iso_seq_without_ptms;
        int sigPSMs;
        double quant;
        Map<String, TranscriptsT> gene_id_map;

        //TODO ||EDITED: Added geneID and mismatches fields||
        // 2021 Expansion
        int allowedMismatches = 0;
        String targetTranscriptID = "";


        String line;
        while ((line = reader.readLine()) != null) {
            // TODO ||Headers exlcuded||
            if ((line.toLowerCase().startsWith("experiment")) || (line.toLowerCase().startsWith("sample"))) {
                continue;
            }
            // TODO ||Tokenize peptide input||
            ArrayList<String> tokens = new ArrayList<>(Arrays.asList(Utils.tokenize(line, "\t", false)));
            //using only the tokens needed.
            // TODO ||Extracts tissue (0) - Label||
            tissue = tokens.get(0).trim();
            // TODO ||Extracts peptide string (1)||
            peptide_string = tokens.get(1).trim();
            // TODO ||Extracts PSMs (2)||
            String sigPsmStr = tokens.get(2).trim();
            if (sigPsmStr.length() == 0) {
                sigPsmStr = "0";
            }
            sigPSMs = Integer.parseInt(sigPsmStr);
            // TODO ||Extracts Quant (3)||
            String quantStr = tokens.get(3).trim();
            if (quantStr.length() == 0) {
                quantStr = "0";
            }
            quant = Double.parseDouble(quantStr);

            // If peptide filter mode is on, and the peptide line has six fields filled.
            if (PepGenomeTool.usePeptideFilter && tokens.size()==6) {
                System.out.println("TabInputPeptidefileParser.read(...) - use filter & six fields");

                // TODO ||EDITED: Extract Allowed Mismatches (4)||
                allowedMismatches = Integer.parseInt(tokens.get(4).trim());


                // TODO ||EDITED: Extract Target Trans ID (5)|| - If left blank, map to all**  Remove testing
                    targetTranscriptID = tokens.get(5).trim(); //working
                    if (targetTranscriptID.equals("")) {
                        targetTranscriptID = "all"; //working
                    }

                //System.out.println(targetTranscriptID); //testing

            } else {
                //TODO ||Remove testing||
                System.out.println("TabInputPeptidefileParser.read(...) - DON'T use filter & six fields");
            }

            //clearing the tokens list for the next iteration.
            tokens.clear();

            // TODO ||Condition: PSMs > 0 - Peptide appears in protein||
            if (sigPSMs > 0) {
                //the matching will only use the amino acids.
                // TODO ||Remove PTMs from the iso sequence, only using amino acids||
                iso_seq_without_ptms = Utils.make_iso_sequence(Utils.remove_ptms(peptide_string));

                // TODO ||Checks coordwrapper if peptide present||
                if (!coordwrapper.isPeptidePresent(iso_seq_without_ptms)) {
                    //the gene_id_map.find_peptide function will match the peptide.
                    // TODO ||Match peptide, produce gene id map using KmerTreeMap/KmerSortedMap - EDITED||
                    // Using peptide filter mode
                    if (PepGenomeTool.usePeptideFilter) {
                        System.out.println("TabInputPeptideFileParser: USING peptide filter --> k.find_peptide(...)");
                        gene_id_map = k.find_peptide(iso_seq_without_ptms, targetTranscriptID, allowedMismatches);
                    }
                    // Default
                    else {
                        System.out.println("TabInputPeptideFileParser: NOT using peptide filter --> k.find_peptide(...)");
                        gene_id_map = k.find_peptide(iso_seq_without_ptms);
                    }

                    //TODO ||Edited: (See k.getIsVariant used below) Extract variant status (check whether peptide contains mismatches)||


                    for (Map.Entry<String, TranscriptsT> it : gene_id_map.entrySet()) {
                        // TODO ||Add peptide to mapping - Creates PeptideEntry object||
                        mapping.add_peptide(coordwrapper, peptide_string, tissue, sigPSMs, gene_id_map.size(), ofs, quant, it, k.getIsVariant());
                    }
                    if (gene_id_map.isEmpty()) {
                        ofs.write(("No-Gene" + "\t" + peptide_string + "\t" + "No-Transcript" + "\t" + "No-genes" + "\t" + tissue + "\t" + sigPSMs + "\t" + quant + "\n").getBytes());
                    }
                } else {
                    //if the peptide already exists its genomic coordinates dont have to be recalculated.
                    //only the tags and PTMs have to be added
                    ArrayList<PeptideEntry> refVec = coordwrapper.get_existing_peptides_at(iso_seq_without_ptms);
                    for (PeptideEntry aRefVec : refVec) {
                        aRefVec.add_peptide(peptide_string, tissue, sigPSMs, quant, k.getIsVariant());
                    }
                }
            }
        }
        ofs.close();
        reader.close();
        ifs.close();
    }

    private void sparkRead(String file, CoordinateWrapper coordwrapper, MappedPeptides mapping, String unmappedoutput, IKmerMap k) throws Exception {

        FileOutputStream ofs = new FileOutputStream(unmappedoutput);

        SparkSession sparkSession = SparkSession.builder().
                master(SparkConfig.getInstance().getMaster())
                .config("spark.ui.enabled", false)
                .appName("pgatk tab input file parser")
                .getOrCreate();

        Dataset<Row> tsv = sparkSession.read()
                .option("sep", "\t")
                .csv(file);

        List<Row> rows = tsv.collectAsList();
        for (Row r : rows) {
            String tissue = r.getString(0).trim();
            if ((tissue.toLowerCase().startsWith("experiment")) || (tissue.toLowerCase().startsWith("sample"))) {
                continue;
            }

            String peptide_string = r.getString(1).trim();
            String sigPsmStr = r.getString(2).trim();
            if (sigPsmStr.length() == 0) {
                sigPsmStr = "0";
            }
            int sigPSMs = Integer.parseInt(sigPsmStr);
            String quantStr = r.getString(3).trim();
            if (quantStr.length() == 0) {
                quantStr = "0";
            }
            double quant = Double.parseDouble(quantStr);

            if (sigPSMs > 0) {
                //the matching will only use the amino acids.
                String iso_seq_without_ptms = Utils.make_iso_sequence(Utils.remove_ptms(peptide_string));

                if (!coordwrapper.isPeptidePresent(iso_seq_without_ptms)) {
                    //the gene_id_map.find_peptide function will match the peptide.
                    Map<String, TranscriptsT> gene_id_map = k.find_peptide(iso_seq_without_ptms);
                    for (Map.Entry<String, TranscriptsT> it : gene_id_map.entrySet()) {
                        mapping.add_peptide(coordwrapper, peptide_string, tissue, sigPSMs, gene_id_map.size(), ofs, quant, it, k.getIsVariant());
                    }
                    if (gene_id_map.isEmpty()) {
                        ofs.write(("No-Gene" + "\t" + peptide_string + "\t" + "No-Transcript" + "\t" + "No-genes" + "\t" + tissue + "\t" + sigPSMs + "\t" + quant + "\n").getBytes());
                    }
                } else {
                    //if the peptide already exists its genomic coordinates dont have to be recalculated.
                    //only the tags and PTMs have to be added
                    ArrayList<PeptideEntry> refVec = coordwrapper.get_existing_peptides_at(iso_seq_without_ptms);
                    for (PeptideEntry aRefVec : refVec) {
                        aRefVec.add_peptide(peptide_string, tissue, sigPSMs, quant, k.getIsVariant());
                    }
                }
            }
        }
        ofs.close();
    }
}