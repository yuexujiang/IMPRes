package kegg;

public class Node {
	private String name;
	private int inDegree=0;
	private int outDegree=0;
	
	Node(String name)
	{
		this.name=name;
	}

	public String getName()
	{
		return this.name;
	}
	
	public int getIndegree()
	{
		return this.inDegree;
	}
	
	public int getOutdegree()
	{
		return this.outDegree;
	}
	
	public void setIndegree(int num)
	{
		this.inDegree=num;
	}
	
	public void setOutdegree(int num)
	{
		this.outDegree=num;
	}
}
