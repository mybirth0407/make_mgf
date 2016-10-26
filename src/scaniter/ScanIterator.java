package scaniter;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;

import modi.Constants;

public abstract class ScanIterator {
  public static int MIN_ASSUMED_CHARGE = 2, MAX_ASSUMED_CHARGE = 3;
  public int sizeOfScans;
  public int scanIndex;
  public String fileName;
  public RandomAccessFile fin;

  public ScanIterator(String fName) throws IOException {
    fileName = fName;
    sizeOfScans = 0;
    scanIndex = 0;
    fin = null;
  }

  public String getFileName() {
    return fileName;
  }

  public int size() {
    return sizeOfScans;
  }

  public int getIndex() {
    return scanIndex;
  }

  public boolean hasNext() throws IOException {
    if (scanIndex < sizeOfScans)
      return true;
    if (fin != null)
      fin.close();
    return false;
  }

  public abstract ArrayList<MSMScan> getNext() throws IOException;

  public static ScanIterator get(String specFile) throws IOException {
    System.out.print("Reading MS/MS spectra.....  ");
    ScanIterator scaniter = null;
    if (specFile.toLowerCase().endsWith(".pkl")) {
      Constants.SPECTRA_FILE_TYPE = Constants.SpectraFileType.PKL;
      scaniter = new PKLIterator(specFile);
    }
    else if (specFile.toLowerCase().endsWith(".mgf")) {
      Constants.SPECTRA_FILE_TYPE = Constants.SpectraFileType.MGF;
      scaniter = new MGFIterator(specFile);
    }
    else if (specFile.toLowerCase().endsWith(".dta")) {
      Constants.SPECTRA_FILE_TYPE = Constants.SpectraFileType.DTA;
      scaniter = new DTAIterator(specFile);
    }
    else if (specFile.toLowerCase().endsWith(".mzxml")) {
      Constants.SPECTRA_FILE_TYPE = Constants.SpectraFileType.MZXML;
      scaniter = new MZXMLIterator(specFile);
    }
    return scaniter;
  }

}
