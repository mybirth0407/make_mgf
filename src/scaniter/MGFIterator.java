package scaniter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class MGFIterator extends ScanIterator {
  public MGFIterator(String fileName) throws IOException {
    super(fileName);

    BufferedReader in = new BufferedReader(new FileReader(fileName));
    String buf;
    int size = 0;
    while ((buf = in.readLine()) != null) {
      if (buf.startsWith("BEGIN IONS"))  // begining of new spectrum
        size++;
    }
    in.close();
    sizeOfScans = size;
    fin = new RandomAccessFile(fileName, "r");
  }

  public ArrayList<MSMScan> getNext() throws IOException {

    ArrayList<MSMScan> scanlist = new ArrayList<MSMScan>();
    MSMScan curScan = null;

    String buf;
    long startOffset = 0, pStart = 0;

    while (true) {
      startOffset = fin.getFilePointer();
      if ((buf = fin.readLine()) == null)
        break;

      // 처음 정보 받는 부분
      if (buf.startsWith("BEGIN")) {  // begining of new spectrum
        String title = String.valueOf(scanIndex + 1);
        int charge = 0;
        double pmz = 0.;
        int scanNo = -1;

        while ((buf = fin.readLine()) != null) {
          if (!buf.contains("=")) {
            break;
          }
          if (buf.startsWith("CHARGE=")) {
            //System.out.println(buf+"ㅋㅋㅋ");
            int st = buf.lastIndexOf('=') + 1;
            int ed = st + 1;
            for (int i = st; i < buf.length(); i++) {
              if (Character.isDigit(buf.charAt(i))) {
                st = i;
                ed = st + 1;
                break;
              }
            }
            for (int i = ed; i < buf.length(); i++) {
              if (!Character.isDigit(buf.charAt(i))) {
                ed = i;
                break;
              }
            }
            charge = Integer.parseInt(buf.substring(st, ed));
          }
          if (buf.startsWith("PEPMASS=")) {
            StringTokenizer
                token =
                new StringTokenizer(buf.substring(buf.indexOf("=") + 1));
            if (token.countTokens() != 0)
              pmz = Double.parseDouble(token.nextToken());
          }
          if (buf.startsWith("SCANS=")) {
            StringTokenizer
                token =
                new StringTokenizer(buf.substring(buf.indexOf("=") + 1));
            if (token.countTokens() != 0)
              scanNo = -1; //only for test.
            //                          scanNo = Integer.parseInt(token.nextToken());
          }
          //추가 08.12
          if (buf.startsWith("TITLE=")) {
            StringTokenizer
                token =
                new StringTokenizer(buf.substring(buf.indexOf("=") + 1));
            if (token.countTokens() != 0)
              title = token.nextToken();
          }
          pStart = fin.getFilePointer(); //mass spectrum start
        }//end while( (buf = fin.readLine()) != null )
        //이제 m/z, intensity받는부분
        if (pmz != 0.) {
          if (charge != 0) {
            if (scanNo != -1)
              curScan = new MSMScan(title, scanNo, pmz, charge);
            else
              curScan = new MSMScan(title, pmz, charge);

            curScan.setOffset(startOffset);
            fin.seek(pStart);
            curScan.readPeakList(fin);
            if (curScan.getSpectrum() != null)
              scanlist.add(curScan);
          }
          else { //charge == 0. 인경우, = 값이 알려지지 않은경우 (제대로 관측되지 않은경우) 2~3까지 charge가 2로도 만들고, 3으로도 만들어서 둘 중 잘되는걸로 맞추려고함
            for (int cs = MIN_ASSUMED_CHARGE; cs <= MAX_ASSUMED_CHARGE; cs++) {
              if (scanNo != -1)
                curScan = new MSMScan(title, scanNo, pmz, cs);//
              else
                curScan = new MSMScan(title, pmz, cs);//

              curScan.setOffset(startOffset);
              fin.seek(pStart);
              curScan.readPeakList(fin);
              if (curScan.getSpectrum() != null)
                scanlist.add(curScan);
            }
          }
        }//end if(pmz!= 0.)
        scanIndex++;
        break; //여기에 break가 있으면 딱 하나의 BEGIN to END ION까지만 읽는것 (메모리때문에?)
      }//end
    }//end while
    return scanlist;
  }
}


