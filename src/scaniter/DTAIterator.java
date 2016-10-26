package scaniter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.StringTokenizer;

import modi.Constants;

public class DTAIterator extends ScanIterator {
  public DTAIterator(String fileName) throws IOException {
    super(fileName);

    BufferedReader in = new BufferedReader(new FileReader(fileName));
    String buf;
    int size = 0;
    boolean newDTA = true;
    while ((buf = in.readLine()) != null) {
      StringTokenizer token = new StringTokenizer(buf);

      if (token.countTokens() != 2) {
        newDTA = true;
        continue;
      }

      if (newDTA && token.countTokens() == 2) {
        size++;
        newDTA = false;
      }
    }
    in.close();

    sizeOfScans = size;
    fin = new RandomAccessFile(fileName, "r");
  }

  public ArrayList<MSMScan> getNext() throws IOException {

    ArrayList<MSMScan> scanlist = new ArrayList<MSMScan>();
    MSMScan curScan = null;

    String buf;
    long startOffset = 0;
    boolean newDTA = true;
    while (true) {
      startOffset = fin.getFilePointer();
      if ((buf = fin.readLine()) == null)
        break;

      StringTokenizer token = new StringTokenizer(buf);
      if (token.countTokens() != 2) {
        newDTA = true;
        continue;
      }

      if (newDTA && token.countTokens() == 2) {
        double precursor = Double.parseDouble(token.nextToken());
        int charge = (int) Double.parseDouble(token.nextToken());
        precursor = (precursor + (charge - 1) * Constants.Proton) / charge;
        curScan = new MSMScan(String.valueOf(++scanIndex), precursor, charge);
        curScan.setOffset(startOffset);
        curScan.readPeakList(fin);
        if (curScan.getSpectrum() != null)
          scanlist.add(curScan);
        break;
      }
    }
    return scanlist;
  }
}
