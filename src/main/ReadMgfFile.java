package main;

import java.io.IOException;
import java.util.ArrayList;

import modi.Constants;
import modi.Spectrum;
import scaniter.MSMScan;
import scaniter.ScanIterator;

/** Author: Yedarm Seong
 * IntelliJ에서는 Wildcard Import시 경로명이 중요함 */
// import org.systemsbiology.jrap.*;
import org.systemsbiology.jrap.stax.*;


public class ReadMgfFile {
  /**
   * Q. 오차 허용?
   */
  static public double tolerance = 0.05;
  /**
   * Q. 허용 범위?
   */
  static public double CorelationRange = 75 * tolerance;
  /**
   * Q. False인데 True일 확률?
   */
  static public double originFT = 0.05;

  public static void main(String[] args) throws IOException {
    /** Q. 펩타이드 조각의 허용 오차?
     * A. MS/MS 의 허용 오차, 장비의 한계 때문에 생김
     * 오차 이내의 범위에서는 구별할 수 없음*/
    Constants.fragmentTolerance = tolerance;
    /** minNumOfPeaksInWindow = 3 일 때 */
    Constants.adjustParametersForInstrument(0);
    Constants.iTRAQSearch = false;

    String testFile = "testFile.mgf";
//    ArrayList<MSMScan> spectrumList = new ArrayList<MSMScan>();
//    spectrumList = ReadMGF(testFile);
    ArrayList<MSMScan> spectrumList = ReadMGF(testFile);
  }

  private static ArrayList<MSMScan> ReadMGF(String filename) throws IOException {
    int readCounter = 0;

    System.out.println("Reading prepxrocessed files," + filename);
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
      /**
       *  Author: Yedarm Seong
       * 더 읽을 게 있다면 반복
       */
      while (scaniter.hasNext()) {
        /** print spectrum count */
        System.out.println(++readCounter);

        ArrayList<MSMScan> scanlist = scaniter.getNext();
        for (MSMScan scan : scanlist) {
          /** Peak Normalization */
          spec = scan.getSpectrum();
          spec.normalizeIntensityLocally();

          /** Peak picking */
          /**
           * Author: Yedarm Seong
           * 2가 이온 이상이고,
           * Q. 인스트럭션 타입이 무엇?
           * A. 장비에 따라 달라질 수 있는 것을 고려하는 것 */
          int extra =
              (spec.getCharge() > 2 && Constants.INSTRUMENTS_TYPE != 0) ? 2 : 0;
          /**
           * Author: Yedarm Seong
           * Q. extra 를 더하는 이유
           * A. 장비 차이에 의해 달라질 수 있는 값들을 보정하기 위해 */
          spec.peakSelection(Constants.selectionWindowSize,
              Constants.minNumOfPeaksInWindow + extra);

          // now we have pre-processed spectrum, spec
          // TODO: write a new preprocessed spectrum file.

          speclist.add(scan);
        }
      }
    } catch (IOException e1) {
      e1.printStackTrace();
    }

    return speclist;
  }
}
