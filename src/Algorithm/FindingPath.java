package Algorithm;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
//import org.apache.commons.math3.stat.descriptive.moment.Variance;



public class FindingPath 
{
	public double dijkstra_findmin(Network net, ArrayList<Node> confidenceSet, ArrayList<Node> startPoint, ArrayList<Node> endPoint, double per, Hashtable<String,Double> ht, int alpha) throws IOException
	{
		System.out.println("startpoint num:"+startPoint.size());
		System.out.println("endPoint num:"+endPoint.size());
		int endnum=endPoint.size();
	
		ArrayList<Node> q =new ArrayList<Node>();
				
			
			for (Node v : net.getNodes().values())
			{
				if(startPoint.contains(v))
				{
					v.setPrevious(null);
					v.setValue(0.0);
					v.setNum(0);
					v.setChildren(new ArrayList<Node>());
				}
				else 
				{
					v.setPrevious(null);
					v.setValue(Double.MAX_VALUE);
					v.setNum(-1);
					v.setChildren(new ArrayList<Node>());
				}       
		        q.add(v);
			}
		
			
			return dijkstra(q,net,endPoint,per,endnum,ht,alpha);
			
			
	}
	
	public double dijkstra(ArrayList<Node> q,Network net, ArrayList<Node> dif,double per,int endnum,Hashtable<String,Double> ht, int alpha)
	{	
		Node u;
		double cover=(double)((endnum-dif.size())/endnum);
		while((!q.isEmpty())&&(cover < per))
		{	
			u=pullLowest(q);
		
			if (u.getValue() == Double.MAX_VALUE) 
			{
				System.out.println("network break");
				break;
			}
			if(dif.contains(u.getBelong()))
			{
				dif.remove(u.getBelong());
	//			System.out.println("FindingPath.java: "+u.getBelong().getname());
				
			}
			for(Node v : u.getDown())
			{
				if(!q.contains(v))
				{
					continue;
				}
				double temp=Double.MAX_VALUE;
				if(alpha==1)
					temp=targetFunc1(u,v,net,ht);
				else if(alpha==2)				
					temp=targetFunc2(u,v,net,ht);
				else if(alpha==3)
					temp=targetFunc3(u,v,net,ht);
				if (temp < v.getValue()) 
				{ 
					q.remove(v);
		            v.setValue(temp);
		            
		            if(v.getPrevious()==null)
		            {
		            	v.setPrevious(u);
		            	u.getChildren().add(v);
		            }
		            else
		            {
		            	Node t=v.getPrevious();
		            	t.getChildren().remove(v);
		            	v.setPrevious(u);
		            	u.getChildren().add(v);
		            }
		            
		            if((u.getBelong()==v)||(v.getBelong()==u))
		            {
		            	v.setNum(u.getNum());
		            }
		            else {
						v.setNum(u.getNum()+1);
					}
	/////////////////////////////////////dream///////////////////////////////////	            
//					if(!v.isIsdif())
//					{
//						v.setExpressionfromV(u.getExpression(1), 1);
//					}
	////////////////////////////////////////////////////////////////////////////	            
		            q.add(v);
		        } 
			}
			
			double newcover=(double)(endnum-dif.size())/endnum;
			cover=newcover;

		}
		System.out.println("cover rate is: "+cover);
		return cover;
	}
	
	public Node pullLowest(ArrayList<Node> q)
	{
		Node result=new Node();
		double value=Double.MAX_VALUE;
		for(Node nd : q)
		{
			if(nd.getValue()<=value)
			{
				result=nd;
				value=nd.getValue();
			}
		}
		q.remove(result);
		return result;
	}
	
	public double targetFunc_unisam(Node S, Node E, Network net, Hashtable<String,Double> ht,int alpha) {
		double value = 0.0;
		String endnodes = S.getname() + E.getname();
		Interaction inter = net.getByendnames(endnodes);
		String pathname=E.getPathway();
		if (inter.getType().equals("type")) {
			value = S.getValue();
		} else {
			if(ht.containsKey(pathname))
			{
				
				if(ht.get(pathname).isInfinite())
				{
					
					value = S.getValue() + E.getWeight()*alpha*69;
				}
				else
				{
					value = S.getValue() + E.getWeight()*alpha*ht.get(pathname);
	//				System.out.println(E.getname()+" "+E.getWeight()+" "+ht.get(pathname));
				}
				
			}
			else
			{
				value = S.getValue() + E.getWeight()*alpha*69;
			}
			
		}

		return value;
	}
//	public double targetFuncdream(Node S, Node E, Network net, Hashtable<String,Double> ht, int alpha)
//	{
//		double value=0.0;
//		String endnodes = S.getname() + E.getname();
//		Interaction inter = net.getByendnames(endnodes);
//		String pathname=E.getPathway();
//		if (inter.getType().equals("type"))
//		{
//			value = S.getValue();
//		}
//		else
//		{
//			double pa=0.0;
//			if(ht.containsKey(pathname))
//			{
//				if(ht.get(pathname).isInfinite())
//				{
//					pa=alpha*69;
//				}
//				else
//				{
//					pa=alpha*ht.get(pathname);
//				}
//			}
//			else
//			{
//				pa=alpha*69;
//			}
//			
//			double[] n1=S.getExpression(1);
//			double[] n2=E.getExpression(1);
//			Variance v=new Variance();
//			if((v.evaluate(n1)==0)||(v.evaluate(n2)==0))
//			{
//				inter.setWeight(1*inter.getWeight());
//			}
//			else
//			{
//				double k=new PearsonsCorrelation().correlation(n1,n2);
//				k=1-Math.abs(k);
//				inter.setWeight(k*inter.getWeight());
//			}
//			
//			value=S.getValue()+1+pa*inter.getWeight();
//		}
//		return value;
//	}
	
	public double targetFunc1(Node S, Node E, Network net, Hashtable<String,Double> ht) {
		double value = 0.0;
		String endnodes = S.getname() + E.getname();
		Interaction inter = net.getByendnames(endnodes);
		String type=inter.getType();
		Pattern p11=Pattern.compile("activation", Pattern.CASE_INSENSITIVE);
		Matcher m11=p11.matcher(type);
		Pattern p12=Pattern.compile("expression", Pattern.CASE_INSENSITIVE);
		Matcher m12=p12.matcher(type);
		Pattern p21=Pattern.compile("inhibition", Pattern.CASE_INSENSITIVE);
		Matcher m21=p21.matcher(type);
		Pattern p22=Pattern.compile("repression", Pattern.CASE_INSENSITIVE);
		Matcher m22=p22.matcher(type);
		
		
		
		
		
		String pathname=E.getPathway();
		if (inter.getType().equals("type")) {
			value = S.getValue();
		} else 
		{
			double nd_penalty=E.getWeight();
			double pathway_penalty=0.0;
			if(ht.containsKey(pathname))
			{
				
				if(ht.get(pathname).isInfinite())
				{
					
					pathway_penalty=69;
				}
				else
				{
					pathway_penalty=ht.get(pathname);
				}
				
			}
			else
			{
				pathway_penalty=69;
			}
			
			if((m11.find())||(m12.find()))
			{
				double fd_product=S.getFc()*E.getFc();
				
				value=1/(1+Math.pow(2.71828, fd_product))*(nd_penalty+pathway_penalty); //simple
			//	value=nd_penalty+pathway_penalty; //more simple
			//	value=pathway_penalty;//path only
			//	value=1/(1+Math.pow(2.71828, fd_product))*(1/(1+Math.pow(2.71828, -1*nd_penalty))+30/(1+Math.pow(2.71828, -1*pathway_penalty)));
				value=value*inter.getWeight();
			}
			else if((m21.find())||(m22.find()))
			{
				double fd_product=S.getFc()*E.getFc();
				value=1/(1+Math.pow(2.71828, -1*fd_product))*(nd_penalty+pathway_penalty);
		//		value=nd_penalty+pathway_penalty;
		//		value=pathway_penalty;
		//		value=1/(1+Math.pow(2.71828, -1*fd_product))*(1/(1+Math.pow(2.71828, -1*nd_penalty))+30/(1+Math.pow(2.71828, -1*pathway_penalty)));
				value=value*inter.getWeight();
			}
			else
			{
				value=0.5*(nd_penalty+pathway_penalty);
		//		value=nd_penalty+pathway_penalty;
		//		value=pathway_penalty;
		//		value=0.5*(1/(1+Math.pow(2.71828, -1*nd_penalty))+30/(1+Math.pow(2.71828, -1*pathway_penalty)));
				value=value*inter.getWeight();
				
			}
			value=S.getValue()+value;

		}
		
		return value;
	}
	public double targetFunc3(Node S, Node E, Network net, Hashtable<String, Double> ht)
	{
		double value = 0.0;
		String endnodes = S.getname() + E.getname();
		Interaction inter = net.getByendnames(endnodes);
		String type=inter.getType();
		if (inter.getType().equals("type")) {
			value = S.getValue();
		} else
		{
			value=S.getValue()+1;
		}
		return value;
	}
	public double targetFunc2new(Node S, Node E, Network net, Hashtable<String,Double> ht)
	{
		double value = 0.0;
		String endnodes = S.getname() + E.getname();
		Interaction inter = net.getByendnames(endnodes);
		
		String pathname=E.getPathway();
		if (inter.getType().equals("type")) {
			value = S.getValue();
		} else {
			double nd_penalty=E.getWeight();
			double pathway_penalty=0.0;
			if(ht.containsKey(pathname))
			{
				
				if(ht.get(pathname).isInfinite())
				{
					
					pathway_penalty=69;
				}
				else
				{
					pathway_penalty=ht.get(pathname);
				}
				
			}
			else
			{
				pathway_penalty=69;
			}
			
			value=S.getValue()+(1/(1+Math.pow(2.71828, -1*inter.getWeight()))+1/(1+Math.pow(2.71828, -1*pathway_penalty)));
			
		}
		
		return value;
	}
	
	public double targetFunc2(Node S, Node E, Network net, Hashtable<String,Double> ht)
	{
		double value = 0.0;
		String endnodes = S.getname() + E.getname();
		Interaction inter = net.getByendnames(endnodes);
		
		String pathname=E.getPathway();
		if (inter.getType().equals("type")) {
			value = S.getValue();
		} else {
			if(ht.containsKey(pathname))
			{
				
				if(ht.get(pathname).isInfinite())
				{
					
					value = S.getValue() + inter.getWeight()*69; //1/(-log2(0.99))
				}
				else
				{
					value = S.getValue() + inter.getWeight()*ht.get(pathname);
				}
				
			}
			else
			{
				value = S.getValue() + inter.getWeight()*69;
			}
		}
		
		return value;
	}
	
	public double targetFunc_step(Node S, Node E, Network net, Hashtable<String,Double> ht,int alpha) {
		double value = 0.0;
		String endnodes = S.getname() + E.getname();
		Interaction inter = net.getByendnames(endnodes);
		String pathname=E.getPathway();
		if (inter.getType().equals("type")) {
			value = S.getValue();
		} else {
			
				value = S.getValue() + inter.getWeight();
			
			
		}

		return value;
	}
	

}
