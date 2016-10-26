package scaniter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;

import org.systemsbiology.jrap.stax.MSXMLParser;
import org.systemsbiology.jrap.stax.Scan;
import org.systemsbiology.jrap.stax.ScanHeader;

public class MZXMLIterator extends ScanIterator {

	private MSXMLParser parser;
	private int maxScanNo, scanNo;
	private Map<Integer, Long> offsetMap;
	
	public MZXMLIterator( String fileName ) throws IOException {
		
		super(fileName);	
		try{
		parser = new MSXMLParser(fileName);	
		maxScanNo = parser.getMaxScanNumber();
	
		int size =0;	
		for (int i=1; i<=maxScanNo; i++){
			Scan scan = parser.rap(i);	
			if( scan == null ) continue;			
 			ScanHeader shead = scan.getHeader();
			if( shead.getMsLevel() > 1 ) size++; // skip MS spectrum(select MSMS spectrum only)						
		}//*/

		offsetMap = parser.getOffsets();	
		sizeOfScans = size;
		scanNo = 1; //necessary
		} catch (Exception e) {
			System.out.println("Abnormal Termination");
			System.exit(1);
		}
	}

	public ArrayList<MSMScan> getNext() throws IOException {	
		
		ArrayList<MSMScan> scanlist = new ArrayList<MSMScan>();
		MSMScan curScan = null;

		for (int i=scanNo; i<=maxScanNo; i++)
		{
			Scan scan = parser.rap(i);	
			if( scan == null ) continue;
			ScanHeader shead = scan.getHeader();
			if( shead.getMsLevel() < 2 ) continue; // skip MS spectrum(select MSMS spectrum only)

			double[][] peakList = scan.getMassIntensityList(); // get peak list array
		
			int precursorCharge = shead.getPrecursorCharge(); // if charge is not determined, guess charge
			if ( precursorCharge < 1 ) {
				// Code below is copied from MsXML2Other of Glossolalia(http://sashimi.sourceforge.net/software_glossolalia.html)
				int intBelow = 0, intAbove = 0;
				for (int j=0; j<peakList[0].length; j++)
					if ( peakList[0][j] <= shead.getPrecursorMz() ) intBelow++; else intAbove++;
				
				if ( (double)intAbove < 0.1 * (double)intBelow ) precursorCharge = 1;	
			}
			
			if( precursorCharge < 1 ) {			
				for(int cs=MIN_ASSUMED_CHARGE; cs<=MAX_ASSUMED_CHARGE; cs++){						
					curScan = new MSMScan(String.valueOf(scanIndex+1), -1, shead.getPrecursorMz(), cs);
					curScan.setOffset( offsetMap.get(i) );			
					curScan.readPeakList( peakList );
					if( curScan.getSpectrum() != null ) scanlist.add(curScan);
				}
			}
			else{			
				curScan = new MSMScan( String.valueOf(scanIndex+1), -1, shead.getPrecursorMz(), precursorCharge );	
				curScan.setOffset( offsetMap.get(i) );	
				curScan.readPeakList( peakList );			
				if( curScan.getSpectrum() != null ) scanlist.add(curScan);
			}
			scanNo = i+1;
			scanIndex++;
			break;
		}
		
		return scanlist;
	}
}
