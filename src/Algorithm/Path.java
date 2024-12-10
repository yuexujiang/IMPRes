package Algorithm;

import java.util.ArrayList;



public class Path {
	private ArrayList<Node> nodes=new ArrayList<Node>();
	private int confinum=0;
	private double value=0.0;
	
	public void setfromarray(String str,Network net)
	{
		String[] strarry=str.split("\t");
		int num=strarry.length;
		for(int i=0;i<num-2;i++)
		{
			Node nd=net.getByName(strarry[i]);
			this.nodes.add(nd);
		}
		this.confinum=Integer.valueOf(strarry[num-2]);
		this.value=Double.valueOf(strarry[num-1]);
	}
	public ArrayList<Node> getNodes()
	{
		return this.nodes;
	}
	public ArrayList<Node> setPath(Network net, Node end, ArrayList<Node> conf)
	{
		if(conf.contains(end.getBelong()))
		{
			this.confinum++;
		}
		if(end==end.getPrevious())
		{
			this.nodes.add(end);
		}
		else if(end.getPrevious()==null)
		{
			this.nodes.add(end);
		}
		else
		{
			this.nodes.add(end);
			this.setPath(net, end.getPrevious(), conf);
			String endnames=end.getPrevious().getname()+end.getname();
			Interaction interab=net.getByendnames(endnames);
			interab.setFreq(interab.getFreq()+1);
		}
		
		return this.nodes;
	}
	public int getConfinum()
	{
		return this.confinum;
	}
	public void setValue(double val)
	{
		this.value=val;
	}
	public double getValue()
	{
		return this.value;
	}
	

}
