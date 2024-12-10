package Algorithm;



public class Interaction 
{
	private Node node1;
	private Node node2;
	private int node1Index;
	private int node2Index;
	private String type; 
	private double weight;
	private Node compound;
	private int freq;
	private double logfreq;
	
	public double getLogfreq() {
		return logfreq;
	}
	public void setLogfreq(double logfreq) {
		this.logfreq = logfreq;
	}
	public int getFreq() {
		return freq;
	}
	public void setFreq(int freq) {
		this.freq = freq;
	}
	Interaction(Node a, Node b)
	{
		this.logfreq=0;
		this.freq=0;
		this.node1=a;
		this.node2=b;
		this.node1Index=a.getIndex();
		this.node2Index=b.getIndex();
		this.type="";
		this.weight=0.0;
		this.compound=new Node();
	}
	public void setCompound(Node nd)
	{
		this.compound=nd;
	}
	public Node getCompound()
	{
		return this.compound;
	}
	public void setWeight(double wt)
	{
		this.weight=wt;
	}
	public double getWeight()
	{
		return this.weight;
	}
	public void setType(String tp)
	{
		this.type=tp;
	}
	public String getType()
	{
		return this.type;
	}
	public Node getNode1()
	{
		return this.node1;
	}
	public Node getNode2()
	{
		return this.node2;
	}
}
