package Algorithm;

import java.util.ArrayList;
import java.util.Vector;

import org.apache.commons.math3.optim.nonlinear.vector.Weight;

public class Node
{
	private double fc;
	private double pvalue;
	private double weight;
	private boolean isdif=false;
	private String fullname;
	private String pathway;	
	private int index;
	private String type; //abst, gene, compound, group
	private double UNIsample;
	private double weightflag;
	private Node previous;
	private ArrayList<Node> children;

	private double value;
	private int num;
	private ArrayList<Node> up;
	private ArrayList<Node> down;
	private ArrayList<Node> contain;
	private Node belong;
	
	private double[] expression1;
	private double[] expression2;
	private ArrayList<Node> truechildren;
	private Node truefather;
	private ArrayList<String> truepathway;
	
//	private double[] expseries;
	
	Node(String name, int num1, int num2)
	{
		this.isdif=false;
		this.pvalue=0.5;
		this.weight=10;
		this.fc=0.0;  //log2
		this.fullname=name;
		this.pathway="";
		this.index=-1;
		this.type="";
//		this.flag="";
		this.children=new ArrayList<Node>();
		this.weightflag=1;
		this.previous=null;
		this.value=0;
		this.num=0;
		this.up=new ArrayList<Node>();
		this.down=new ArrayList<Node>();
		this.contain=new ArrayList<Node>();
		this.belong=null;
		this.expression1=new double[num1];
		this.expression2=new double[num2];
		this.truechildren=new ArrayList<Node>();
		this.truefather=null;
		this.truepathway=new ArrayList<String>();
//		this.expseries=new double[num];
	}
	Node() 
	{
		this.isdif=false;
		this.pvalue=0.5;
		this.fc=0.0;
		this.weight=10;
		this.fullname="";
		this.pathway="";
		this.index=-1;
		this.type="";
//		this.flag="";
		this.children=new ArrayList<Node>();
		this.weightflag=1;
		this.previous=null;
		this.value=0;
		this.up=new ArrayList<Node>();
		this.down=new ArrayList<Node>();
		this.contain=new ArrayList<Node>();
		this.num=0;
		this.belong=null;
		this.truechildren=new ArrayList<Node>();
		this.truefather=null;
		this.truepathway=new ArrayList<String>();
//		this.expression1=new double[num1];
//		this.expression2=new double[num2];
	}
	public void deepcopyclass(int k,Node source_class, Network net)
	{
	//	Node target_class=new Node();
		this.assign(k,source_class);
	//	children
		
	//	previous
	//	value
		for(Node source_inst : source_class.getContain())
		{
			Node target_inst=new Node("any",this.getExpression(1).length,this.getExpression(2).length);		
			target_inst.deepcopyinst(k, source_inst, net);
			
			this.addContain(target_inst);
			target_inst.setBelong(this);
			this.addUP(target_inst);
			this.addDown(target_inst);
			target_inst.addUP(this);
			target_inst.addDown(this);
			Interaction interab=new Interaction(this,target_inst);
			interab.setWeight(0.0);
			interab.setType("type");
			net.addInteraction(interab);
			interab=new Interaction(target_inst,this);
			interab.setWeight(0.0);
			interab.setType("type");
			net.addInteraction(interab);
			
			
			net.addNode(target_inst);
		}
		net.addNode(this);
		
	}
	public void deepcopyinst(int k, Node source_inst, Network net)
	{
	//	Node target_inst=new Node();
		this.assign(k, source_inst);
		for(Node any : source_inst.getUp())
		{
			if(!any.getType().equals("abst"))
			{
				this.addUP(any);
				any.addDown(this);
				
				String type=net.getByendnames(any.getname()+source_inst.getname()).getType();
				Interaction interab=new Interaction(any,this);
				interab.setWeight(1.0);
				interab.setType(type);
				net.addInteraction(interab);
			}
		}
		for(Node any : source_inst.getDown())
		{
			if(!any.getType().equals("abst"))
			{
				this.addDown(any);
				any.addUP(this);
				
				
				String type=net.getByendnames(source_inst.getname()+any.getname()).getType();
				Interaction interab=new Interaction(this,any);
				interab.setWeight(1.0);
				interab.setType(type);
				net.addInteraction(interab);
			}
		}
		
	}
	public void assign(int k, Node n)
	{
		this.setIsdif(n.isIsdif());
		this.setPvalue(n.getPvalue());
		this.setFc(n.getFc());
		this.setWeight(n.getWeight());
		this.setname(n.getname()+"?"+k);
		this.setPathway(n.getPathway());
		this.setIndex(n.getIndex());
		this.setType(n.getType());
		this.setWeightflag(n.getWeightflag());
	}
	
	public ArrayList<String> getTruePathway() {
		return truepathway;
	}
	public void addTruePathway(String a) {
		this.truepathway.add(a);
	}
	
	public Node getTruefather() {
		return truefather;
	}
	public void setTruefather(Node a) {
		this.truefather=a;
	}
	
	public ArrayList<Node> getTrueChildren() {
		return truechildren;
	}
	public void addTrueChildren(Node a) {
		this.truechildren.add(a);
	}
	
	public ArrayList<Node> getChildren() {
		return children;
	}
	public void setChildren(ArrayList<Node> children) {
		this.children = children;
	}
	public boolean isIsdif() {
		return isdif;
	}
	public void setIsdif(boolean isdif) {
		this.isdif = isdif;
	}
	public double getWeight() {
		return weight;
	}
	public void setWeight(double weight) {
		this.weight = weight;
	}
	public double getPvalue() {
		return pvalue;
	}
	public void setPvalue(double pvalue) {
		this.pvalue = pvalue;
	}
	public double getFc() {
		return fc;
	}
	public void setFc(double fc) {
		this.fc = fc;
	}
	
	
	public double[] getExpression(int i)
	{
		double[] expression={};
		if(i==1)
		{
			expression=this.expression1;
		}
		if(i==2)
		{
			expression=this.expression2;
		}
//		if(i==3)
//		{
//			expression=this.expseries;
//		}
		return expression;
	}
	public void setExpressionfromarray(String[] db, int i)
	{
		int s=db.length;
		if(i==1)
		{
			for(int j=0;j<s;j++)
			{
				this.expression1[j]=Double.valueOf(db[j]);
			}
		}
		if(i==2)
		{
			for(int j=0;j<s;j++)
			{
		//		System.out.println(db[j]);
				this.expression2[j]=Double.valueOf(db[j]);
			}
		}
//		if(i==3)
//		{
//			for(int j=0;j<s;j++)
//			{
//				this.expseries[j]=Double.valueOf(db[j]);
//			}
//		}
	}
	public void setExpressionfromV(double[] db, int i)
	{
		if(i==1)
		{
			this.expression1=db;
		}
		if(i==2)
		{
			this.expression2=db;
		}
//		if(i==3)
//		{
//			this.expseries=db;
//		}
	}
	public void setExpression4group(int i, int num, Network net)
	{	
		String[] str_out=this.fullname.split("\\|com\\|");
		double[] vec_out=new double[num];
		for(int j=0;j<num;j++)
		{
			vec_out[j]=Double.MAX_VALUE;
		}
		for(int j=1;j<str_out.length;j++)
		{
			String[] str_in=str_out[j].split("_");
			double[] vec_in=new double[num];
			for(int k=0;k<num;k++)
			{
				double db=0.0;
				for(int z=0;z<str_in.length;z++)
				{
					Node in=net.getByName("abst#"+str_in[z]);
			//		if(in.getExpression(i).size()==0)
			//		{
			//			return;
			//		}
					if(in.getExpression(i)[k]>db)
					{
						db=in.getExpression(i)[k];
					}
				}
				vec_in[k]=db;
			}
			for(int p=0;p<num;p++)
			{
				if(vec_in[p]<vec_out[p])
				{
					vec_out[p]=vec_in[p];
					
				}
			}
		}
		if(i==1)
		{
			this.expression1=vec_out;
		}
		if(i==2)
		{
			this.expression2=vec_out;
		}	
//		if(i==3)
//		{
//			this.expseries=vec_out;
//		}
	}
	public void setNum(int num)
	{
		this.num=num;
	}
	public int getNum()
	{
		return this.num;
	}
	public void setValue(double vl)
	{
		this.value=vl;
	}
	public double getValue()
	{
		return this.value;
	}
	public void setPrevious(Node nd)
	{
		this.previous=nd;
	}
	public Node getPrevious()
	{
		return this.previous;
	}
	public void setWeightflag(double wt)
	{
		this.weightflag=wt;
	}
	public double getWeightflag()
	{
		return this.weightflag;
	}
	public void setUnisample(double str)
	{
		this.UNIsample=str;
	}
	public double getunisample(double str)
	{
		return this.UNIsample;
	}
	public void addUP(Node nd)
	{
		this.up.add(nd);
	}
	public void addDown(Node nd)
	{
		this.down.add(nd);
	}
	public ArrayList<Node> getDown()
	{
		return this.down;
	}
	public ArrayList<Node> getUp()
	{
		return this.up;
	}
	public void addContain(Node nd)
	{
		this.contain.add(nd);
	}
	public ArrayList<Node> getContain()
	{
		return this.contain;
	}
	public void setBelong(Node nd)
	{
		this.belong=nd;
	}
	public Node getBelong()
	{
		return this.belong;
	}
	public void setname(String name)
	{
		this.fullname=name;
	}
	public String getname()
	{
		return this.fullname;
	}
	public void setPathway(String pname)
	{
		this.pathway=pname;
	}
	public String getPathway()
	{
		return this.pathway;
	}
	public void setIndex(int num)
	{
		this.index=num;
	}
	public int getIndex()
	{
		return this.index;
	}
	public void setType(String tp)
	{
		this.type=tp;
	}
	public String getType()
	{
		return this.type;
	}
	
}
