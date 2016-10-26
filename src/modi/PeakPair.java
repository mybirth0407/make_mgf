package modi;

public class PeakPair {
	public Peak first;
	public Peak second;
	
	public PeakPair(Peak first, Peak second)
	{
		this.first = first;
		this.second = second;
	}
	public	Peak	getFirst()	{ return first; }
	public	Peak	getSecond()	{ return second; }
	
	public	int		getAttribute()
	{
		if(first != null)
		{
			if(second != null)
				return 0;		// shared
			else
				return 1;		// firstOnly
		}
		else
		{
			if(second != null)
				return 2;		// secondOnly
			else
				return -1;		// no peak : error
		}
	}

	public	String	toString()
	{
		return new String("(" + first + "," + second + ")");
	}
}
