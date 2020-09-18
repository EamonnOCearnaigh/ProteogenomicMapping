package org.bigbio.pgatk.pepgenome.io;

import org.bigbio.pgatk.pepgenome.common.FastaEntry;
import org.bigbio.pgatk.pepgenome.common.Utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;

/**
 * This fastaparser reads fasta files. it will convert input sequences into
 * iso-sequences [I, L] will be converted to J
 *
 * @author ypriverol
 *
 */
public class FastaParser implements Serializable {

    private static final long serialVersionUID = -4871649120608342251L;

    //inputstream.
    private BufferedReader br = null;

    //current line.
    private String mLine = "";

    //meyers singleton instance.
    private static FastaParser instance = null;

    private FastaParser() {
    }

    //singleton get_instance method.
    public static FastaParser get_instance() {
        if (instance == null) {
            instance = new FastaParser();
        }
        return instance;
    }

    //opens the file
    public boolean open(String file) throws Exception {
        if (br == null) {
            br = new BufferedReader(new FileReader(file));
            mLine = "";
        }
        return br.ready();
    }

    //closes the filestream
    public final void close() {
        mLine = "";
        try {
            br.close();
        } catch (IOException e) {
            //e.printStackTrace();
        }
        br = null;
    }

/*
    // TODO Edit FastaParser.nextEntry() - Attached new method, see below. (Working so far)
    //parses and returns the next FASTA entry.
    public final FastaEntry nextEntry() throws IOException {
        if (mLine.isEmpty()) {
            mLine = br.readLine();
        }

        if (mLine == null) {
            return new FastaEntry("", "");
        }
        String header = mLine;
        StringBuilder sequenceBuilder = new StringBuilder();
        while ((mLine = br.readLine()) != null && !mLine.startsWith(">")) {
            if (mLine.startsWith(">")) {
                //extractOffset(mLine);
            }
            else {
                sequenceBuilder.append(Utils.make_iso_sequence(mLine));
            }

        }
        if (mLine == null || !mLine.startsWith(">")) {
            mLine = "";
        }

        return new FastaEntry(header, sequenceBuilder.toString());
    }

 */

    // Copy of original method.
    //parses and returns the next FASTA entry.
    public final FastaEntry nextEntry() throws IOException {
        if (mLine.isEmpty()) {
            mLine = br.readLine();
        }

        if (mLine == null) {
            return new FastaEntry("", "");
        }
        String header = mLine;
        StringBuilder sequenceBuilder = new StringBuilder();
        while ((mLine = br.readLine()) != null && !mLine.startsWith(">")) {
            sequenceBuilder.append(Utils.make_iso_sequence(mLine));
        }
        if (mLine == null || !mLine.startsWith(">")) {
            mLine = "";
        }

        return new FastaEntry(header, sequenceBuilder.toString());
    }


/*
    // TODO Edit - New Method
    private void extractOffset (String line) {
        if (line.contains("spos")) {
            System.out.println(true);
        }
    }

 */

}