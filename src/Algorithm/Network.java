package Algorithm;

import java.util.*;


public class Network {

	public Hashtable<String, Node> nodes;
	private Hashtable<String, Interaction> interactions=new Hashtable<String, Interaction>();
	
	private int exp_num=0;
	
	public int getExp_num() {
		return exp_num;
	}

	public void setExp_num(int exp_num) {
		this.exp_num = exp_num;
	}

	public Network() {
		this.exp_num=0;
		this.nodes=new Hashtable<String, Node>();
	}
	
	public Hashtable<String, Node> getNodes()
	{
		return this.nodes;
	}

	public Node getByName(String nodename)
	{
		return this.nodes.get(nodename);
	}
	public void addNode(Node nd)
	{
		this.nodes.put(nd.getname(), nd);
	}
	public void addInteraction(Interaction edge)
	{
		String node1=edge.getNode1().getname();
		String node2=edge.getNode2().getname();
		String endnodes=node1+node2;
		this.interactions.put(endnodes, edge);
	}
	public Interaction getByendnames(String endnames)
	{
		return this.interactions.get(endnames);
	}
	public Hashtable<String, Interaction> getInteractions()
	{
		return this.interactions;
	}
	public Boolean nodeExist(String nodename)
	{
		if(this.nodes.containsKey(nodename))
		{
			return true;
		}
		return false;
	}
	public Boolean edgeExist(String endnames)
	{
		if(this.interactions.containsKey(endnames))
		{
			return true;
		}
		return false;
	}
}
