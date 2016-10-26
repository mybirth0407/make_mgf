package main;   
import java.io.IOException;
import java.util.ArrayList;

import modi.Constants;
import modi.Spectrum;
import scaniter.MSMScan;
import scaniter.ScanIterator;
import org.systemsbiology.*;


public class ReadMgfFile {

  static double tolerance = 0.05;
  static public double CorelationRange = 75 * tolerance;
  static public double originFT = 0.05;

  public static void main(String[] args) throws IOException {
    Constants.fragmentTolerance = tolerance;
    Constants.adjustParametersForInstrument(0);
    Constants.iTRAQSearch = false;
    
    
    String testFile = "testFile.mgf";
    ArrayList<MSMScan> spectrumList = new ArrayList<MSMScan>();
    spectrumList = ReadMGF(testFile);
    
    
  }

  private static ArrayList<MSMScan> ReadMGF(String filename) throws IOException {
    int readCounter = 0;

    System.out.println("Reading preprocessed files," + filename);
    ScanIterator scaniter = null;
    ArrayList<MSMScan> speclist = new ArrayList<MSMScan>();

    scaniter = ScanIterator.get(filename);

    if (scaniter == null || scaniter.size() == 0) {
      System.out.println("Failed to read msms spectra file");
    }
    // print the number of spectra
    System.out.println(scaniter.size() + " scans");
    
    Spectrum spec;
    try {
      while (scaniter.hasNext()) {
        // print spectrum count
        System.out.println(++readCounter);

        ArrayList<MSMScan> scanlist = scaniter.getNext();
        for (MSMScan scan : scanlist) {
          //Peak Normalization
          spec = scan.getSpectrum();
          spec.normalizeIntensityLocally();
          
          //Peak picking
          int extra = (spec.getCharge() > 2 && Constants.INSTRUMENTS_TYPE != 0) ? 2 : 0;
          spec.peakSelection(Constants. selectionWindowSize, Constants.minNumOfPeaksInWindow + extra);
          
          //now we have pre-processed spectrum, spec
          //TODO: write a new preprocessed spectrum file.
          
          speclist.add(scan);
          
        }
      }
    } catch (IOException e1) {
      e1.printStackTrace();
    }

    return speclist;
  }
}
