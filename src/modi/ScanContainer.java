package modi;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

public class ScanContainer extends ArrayList<ScanCap> {
	
	private int MAXCHARGE = 3;
	public ScanContainer( int cs ){ MAXCHARGE = cs; }
	
	public int parseFromMGF( String fileName ) throws IOException
	{
		if( fileName == null ) {
			System.out.println( "The MS/MS dataset was not specified." );
			return 0;
		}

		RandomAccessFile in= new RandomAccessFile(fileName, "r");
		String buf;
		long startPeakOffset = 0;
		ScanCap curScan = null;
		
		while( (buf = in.readLine()) != null ){
		
			if( buf.startsWith("BEGIN") ){	// begining of new spectrum	
				String title= String.valueOf( this.size()+1 );
				int charge= 0;
				double pmz= 0.;

				while( (buf = in.readLine()) != null )
				{	
					if( !buf.contains("=") ) break;
					
					if( buf.startsWith("TITLE=") ){
						title = buf.substring(buf.indexOf("=")+1).trim();	
						title = title.replaceAll("[/,:*?\"<>|\\\\]", "_");
					}
					if( buf.startsWith("CHARGE=") ){												
						int st = buf.lastIndexOf('=')+1;
						int ed = st+1;
						for(int i=st; i<buf.length(); i++){
							if( Character.isDigit( buf.charAt(i) ) ){
								st = i;
								ed = st+1;
								break;
							}
						}
						for(int i=ed; i<buf.length(); i++){
							if( !Character.isDigit( buf.charAt(i) ) ){
								ed = i;
								break;
							}
						}
						charge= Integer.parseInt( buf.substring(st,ed) );
					}
					if( buf.startsWith("PEPMASS=") ){						
						StringTokenizer token = new StringTokenizer(buf.substring(buf.indexOf("=")+1));
						if( token.countTokens() != 0 )
							pmz = Double.parseDouble( token.nextToken() );
					}
					startPeakOffset = in.getFilePointer();
				}
				
				if( pmz != 0. ){
					if( charge != 0 ){			
						curScan = new ScanCap(title, pmz, charge);				
						curScan.setOffset( startPeakOffset );			
						if( curScan != null ){ this.add(curScan); }
					}
					else{
						for(int cs=2; cs<=MAXCHARGE; cs++){
							curScan = new ScanCap(title+"."+cs, pmz, cs);				
							curScan.setOffset( startPeakOffset );			
							if( curScan != null ){ this.add(curScan); }
						}
					}
				}
			}	
		}
		in.close();
		return 0;
	}

	public int parseFromPKL( String fileName ) throws IOException
	{
		if( fileName == null ) {
			System.out.println( "The MS/MS dataset was not specified." );
			return 0;
		}
				
		RandomAccessFile in= new RandomAccessFile(fileName, "r");
		
		String buf;
		ScanCap curScan = null;
		int scanSize=1;
		while( (buf = in.readLine()) != null ) 
		{
			StringTokenizer token = new StringTokenizer(buf);
			if( token.countTokens() == 3 ){
				
				double precursor = Double.parseDouble(token.nextToken());
				double precursorIntensity = Double.parseDouble(token.nextToken());				
				int charge = (int)Double.parseDouble(token.nextToken());
				
				curScan = new ScanCap( String.format("%d.%d.%d", scanSize, scanSize, charge), precursor, charge );				
				curScan.setOffset( in.getFilePointer() );	
				if( curScan != null ){ this.add(curScan); scanSize++;}
			}
		}	
		in.close();
		return 0;
	}
	
	public int parseFromDTA( String fileName ) throws IOException
	{
		if( fileName == null ) {
			System.out.println( "The MS/MS dataset was not specified." );
			return 0;
		}
		RandomAccessFile in= new RandomAccessFile(fileName, "r");
		
		String buf;
		ScanCap curScan = null;
		boolean newDTA = true;
		int scanSize = 1;
		while( (buf = in.readLine()) != null ) 
		{
			StringTokenizer token = new StringTokenizer(buf);
			
			if( token.countTokens() != 2 ){
				newDTA = true;
				continue;
			}
			
			if( newDTA && token.countTokens() == 2 ){
				double precursor = Double.parseDouble(token.nextToken());
				int charge = (int)Double.parseDouble(token.nextToken());
				precursor = (precursor+(charge-1)*Constants.Proton)/charge;				
				curScan = new ScanCap( String.format("%d.%d.%d", scanSize, scanSize, charge), precursor, charge );				
				curScan.setOffset( in.getFilePointer() );	
				if( curScan != null ){ 
					this.add(curScan); 
					scanSize++;
					newDTA = false; 
				}
			}
		}	
		in.close();
		return 0;
	}
	
	public int parseFromMS2( String fileName ) throws IOException
	{
		if( fileName == null ) {
			System.out.println( "The MS/MS dataset was not specified." );
			return 0;
		}
		String baseName = fileName.replace("\\", "/");
		baseName = baseName.substring(baseName.lastIndexOf('/')+1, baseName.lastIndexOf('.'));
		RandomAccessFile in= new RandomAccessFile(fileName, "r");
		String buf;
		long startPeakOffset = 0;
		ScanCap curScan = null;
		
		while( (buf = in.readLine()) != null ){
			
			if( buf.startsWith("S") ){	// begining of new spectrum	
				String title= String.valueOf( this.size()+1 );				
				StringTokenizer token = new StringTokenizer(buf);
				if( token.countTokens() != 4 ) continue;
				token.nextToken();//S
				String startScan = token.nextToken();
				String endScan = token.nextToken();
				title = baseName+"."+startScan+"."+endScan+".";
				
				double pmz = Double.parseDouble( token.nextToken() );
				int charge= 0;				
				
				while( (buf = in.readLine()) != null ) {	
					if( !Character.isLetter(buf.charAt(0)) ) break;
					if( buf.startsWith("Z") ){
						token = new StringTokenizer(buf);
						if( token.countTokens() != 3 ) break;
						token.nextToken();//Z
						charge = Integer.parseInt( token.nextToken() );
					}
					startPeakOffset = in.getFilePointer();
				}
				
				if( charge != 0 ){
					title += charge;
					curScan = new ScanCap(title, pmz, charge);				
					curScan.setOffset( startPeakOffset );			
					if( curScan != null ){ this.add(curScan); }
				}
			}	
		}
		in.close();
		return 0;
	}
	
	public void sortByMW(){
		if( this.size() < 2 ) return;
		Collections.sort(this);
	}
	
}
