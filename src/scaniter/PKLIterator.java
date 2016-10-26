package scaniter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class PKLIterator  extends ScanIterator {
	public PKLIterator( String fileName ) throws IOException {
		
		super(fileName);
		
		BufferedReader in = new BufferedReader( new FileReader(fileName) );
		String buf;
		int size = 0;			
		while( (buf = in.readLine()) != null ){
			StringTokenizer token = new StringTokenizer(buf);
			if( token.countTokens()== 3 )	// begining of new spectrum
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
		long startOffset = 0;
		while( true ){
			
			startOffset = fin.getFilePointer();
			if( (buf = fin.readLine()) == null ) break;
			StringTokenizer token = new StringTokenizer(buf);
			if( token.countTokens() == 3 ){			
				double precursor = Double.parseDouble(token.nextToken());
				double precursorIntensity = Double.parseDouble(token.nextToken());				
				int charge = (int)Double.parseDouble(token.nextToken());
				
				curScan = new MSMScan( String.valueOf(++scanIndex), precursor, charge );				
				curScan.setOffset( startOffset );
				curScan.readPeakList(fin);
				if( curScan.getSpectrum() != null ) scanlist.add(curScan);
				break;
			}
		}	

		return scanlist;
	}
}
