package Algorithm;


import java.io.*;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.correlation.*;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.descriptive.moment.*;
import org.rosuda.REngine.*;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;

public class IO 
{
	public double trans(double value, double base)
	{
		double v=Math.log(value);
		double b=Math.log(base);
		
		return v/b;
	}
	
	public ArrayList<String> getexistingEdgelist(Network net)
	{
		ArrayList<String> al=new ArrayList<String>();
		for(Map.Entry<String, Interaction> i : net.getInteractions().entrySet())
		{
			Interaction edge=i.getValue();
			if(edge.getType().equals("type"))
			{}
			else
			{
				Node nd1=edge.getNode1();
				Node nd2=edge.getNode2();
				String abst_abst=nd1.getBelong().getname()+"_"+nd2.getBelong().getname();
				if(al.contains(abst_abst))
				{}
				else
				{
					al.add(abst_abst);
				}
				
			}
		}
		
		
		return al;
	}
	
	public Network readTFmRNA4Inte(ArrayList<String> al, Network net,  String Stringname, int num1, int num2) throws IOException
	{
		int jishu=0;
		File StringNET=new File(Stringname);
		BufferedReader reader = new BufferedReader(new FileReader(StringNET));
		String line="";
		while ((line = reader.readLine()) != null)
		{
//			System.out.println("IO: "+line);
			String[] str=line.split("\t");
//			System.out.println(line);
//			if(Integer.valueOf(str[2])<900)
//				continue;
			String forcheck1=str[0];
			String forcheck2=str[2];
			String relation=str[1];
			
			
			if(al.contains(forcheck1+"_"+forcheck2))
			{}
			else if(al.contains(forcheck2+"_"+forcheck1))
			{}
			else
			{
				jishu++;
				String instName1="inst#gene#unknown#"+forcheck1.substring(5);
				if(net.nodeExist(instName1))
				{
					
				}
				else
				{
					if(net.nodeExist(forcheck1))
					{
						Node node1=new Node(instName1,num1,num2);
						node1.setType("gene");
						node1.setPathway("unknown");
						Node abst1=net.getByName(forcheck1);
						net.addNode(node1);
						node1.setBelong(abst1);
						abst1.addContain(node1);
						node1.addDown(abst1);
						node1.addUP(abst1);
						abst1.addDown(node1);
						abst1.addUP(node1);
						Interaction interab=new Interaction(node1,abst1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst1, node1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
					else
					{
						Node node1=new Node(instName1,num1,num2);
						node1.setType("gene");
						node1.setPathway("unknown");
						Node abst1=new Node(forcheck1,num1,num2);
						abst1.setType("abst");
						net.addNode(abst1);
						net.addNode(node1);
						node1.setBelong(abst1);
						abst1.addContain(node1);
						node1.addDown(abst1);
						node1.addUP(abst1);
						abst1.addDown(node1);
						abst1.addUP(node1);
						Interaction interab=new Interaction(node1,abst1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst1, node1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
				}
				String instName2="inst#gene#unknown#"+forcheck2.substring(5);
				if(net.nodeExist(instName2))
				{
					
				}
				else
				{
					if(net.nodeExist(forcheck2))
					{
						Node node2=new Node(instName2,num1,num2);
						node2.setType("gene");
						node2.setPathway("unknown");
						Node abst2=net.getByName(forcheck2);
						net.addNode(node2);
						node2.setBelong(abst2);
						abst2.addContain(node2);
						node2.addDown(abst2);
						node2.addUP(abst2);
						abst2.addDown(node2);
						abst2.addUP(node2);
						Interaction interab=new Interaction(node2,abst2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst2, node2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
					else
					{
						Node node2=new Node(instName2,num1,num2);
						node2.setType("gene");
						node2.setPathway("unknown");
						Node abst2=new Node(forcheck2,num1,num2);
						abst2.setType("abst");
						net.addNode(abst2);
						net.addNode(node2);
						node2.setBelong(abst2);
						abst2.addContain(node2);
						node2.addDown(abst2);
						node2.addUP(abst2);
						abst2.addDown(node2);
						abst2.addUP(node2);
						Interaction interab=new Interaction(node2,abst2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst2, node2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
				}
				Node nd1=net.getByName(instName1);
				Node nd2=net.getByName(instName2);
				Interaction interab=new Interaction(nd1,nd2);
				interab.setType(relation.toLowerCase());
				net.addInteraction(interab);
				interab.setWeight(1.0);
				
				
				nd1.addDown(nd2);
				nd2.addUP(nd1);
				
				al.add(forcheck1+"_"+forcheck2);
				
			}
			
		}
		System.out.println("# of PPI after TF_mRNA: "+jishu);
		System.out.println("# of node after TF_mRNA: "+net.nodes.size());
		return net;
	}
	public Network readPPI4Inte(ArrayList<String> al, Network net,  String Stringname, int num1, int num2) throws IOException
	{
		int jishu=0;
		File StringNET=new File(Stringname);
		BufferedReader reader = new BufferedReader(new FileReader(StringNET));
		String line="";
		while ((line = reader.readLine()) != null)
		{
//			System.out.println("IO: "+line);
			String[] str=line.split(" ");
//			System.out.println(line);
//			if(Integer.valueOf(str[2])<900)
//				continue;
			String forcheck1="abst#"+str[0];
			String forcheck2="abst#"+str[1];
			
			
			if(al.contains(forcheck1+"_"+forcheck2))
			{}
			else if(al.contains(forcheck2+"_"+forcheck1))
			{}
			else
			{
				jishu++;
				String instName1="inst#gene#unknown#"+forcheck1.substring(5);
				if(net.nodeExist(instName1))
				{
					
				}
				else
				{
					if(net.nodeExist(forcheck1))
					{
						Node node1=new Node(instName1,num1,num2);
						node1.setType("gene");
						node1.setPathway("unknown");
						Node abst1=net.getByName(forcheck1);
						net.addNode(node1);
						node1.setBelong(abst1);
						abst1.addContain(node1);
						node1.addDown(abst1);
						node1.addUP(abst1);
						abst1.addDown(node1);
						abst1.addUP(node1);
						Interaction interab=new Interaction(node1,abst1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst1, node1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
					else
					{
						Node node1=new Node(instName1,num1,num2);
						node1.setType("gene");
						node1.setPathway("unknown");
						Node abst1=new Node(forcheck1,num1,num2);
						abst1.setType("abst");
						net.addNode(abst1);
						net.addNode(node1);
						node1.setBelong(abst1);
						abst1.addContain(node1);
						node1.addDown(abst1);
						node1.addUP(abst1);
						abst1.addDown(node1);
						abst1.addUP(node1);
						Interaction interab=new Interaction(node1,abst1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst1, node1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
				}
				String instName2="inst#gene#unknown#"+forcheck2.substring(5);
				if(net.nodeExist(instName2))
				{
					
				}
				else
				{
					if(net.nodeExist(forcheck2))
					{
						Node node2=new Node(instName2,num1,num2);
						node2.setType("gene");
						node2.setPathway("unknown");
						Node abst2=net.getByName(forcheck2);
						net.addNode(node2);
						node2.setBelong(abst2);
						abst2.addContain(node2);
						node2.addDown(abst2);
						node2.addUP(abst2);
						abst2.addDown(node2);
						abst2.addUP(node2);
						Interaction interab=new Interaction(node2,abst2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst2, node2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
					else
					{
						Node node2=new Node(instName2,num1,num2);
						node2.setType("gene");
						node2.setPathway("unknown");
						Node abst2=new Node(forcheck2,num1,num2);
						abst2.setType("abst");
						net.addNode(abst2);
						net.addNode(node2);
						node2.setBelong(abst2);
						abst2.addContain(node2);
						node2.addDown(abst2);
						node2.addUP(abst2);
						abst2.addDown(node2);
						abst2.addUP(node2);
						Interaction interab=new Interaction(node2,abst2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst2, node2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
				}
				Node nd1=net.getByName(instName1);
				Node nd2=net.getByName(instName2);
				Interaction interab=new Interaction(nd1,nd2);
				interab.setType("PPI");
				net.addInteraction(interab);
				interab.setWeight(1.0);
				
				interab=new Interaction(nd2, nd1);
				interab.setType("PPI");
				net.addInteraction(interab);
				interab.setWeight(1.0);
				
				nd1.addDown(nd2);
				nd1.addUP(nd2);
				nd2.addDown(nd1);
				nd2.addUP(nd1);
				
				al.add(forcheck1+"_"+forcheck2);
				
			}
			
		}
		System.out.println("# of PPI : "+jishu);
		System.out.println("# of node : "+net.nodes.size());
		return net;
	}
	
	public Network readPPIfromFile(ArrayList<String> al, Network net,  String Stringname, String Stringmap, int num1, int num2) throws IOException
	{
		Hashtable<String, String> ht_dipmapABST=new Hashtable<String, String>();
	//	String forcheck="";
		int jishu=0;
		
//		File dipmapfile=new File(DIPmap);
//		BufferedReader reader = new BufferedReader(new FileReader(dipmapfile));
//		String line="";
//		while ((line = reader.readLine()) != null)
//		{
//			String[] str=line.split("\t");
//			if(str.length<2)
//			{continue;}
//			else
//			{
//				String dipID=str[0];
//				String keggid="abst#"+str[1].split(";")[0];   //seperated by ";"
//				ht_dipmapABST.put(dipID, keggid);
//			}
//			
//		}
//		
//		File dipNET=new File(DIPname);
//		reader = new BufferedReader(new FileReader(dipNET));
//		while ((line = reader.readLine()) != null)
//		{
//		//	System.out.println(line);
//			String[] str=line.split("\t");
//			String forcheck1="";
//			String forcheck2="";
//			Pattern p1=Pattern.compile("(DIP-\\d*N)", Pattern.CASE_INSENSITIVE);
//			Matcher m1=p1.matcher(str[0]);
//			int duo=0;
//			if(m1.find())
//			{
//				String id1=m1.group(1);
//				if(ht_dipmapABST.containsKey(id1))
//				{
//					forcheck1=ht_dipmapABST.get(id1);
//				}
//				else
//				{
//					forcheck1="abst#"+id1;
//				}
//				duo++;
//			}
//			Matcher m2=p1.matcher(str[1]);
//			
//			if(m2.find())
//			{
//				String id2=m2.group(1);
//				if(ht_dipmapABST.containsKey(id2))
//				{
//					forcheck2=ht_dipmapABST.get(id2);
//				}
//				else
//				{
//					forcheck2="abst#"+id2;
//				}
//				duo++;
//			}
//			
//			if(duo!=2)
//			{
//				continue;
//			}
//			
//			if(al.contains(forcheck1+"_"+forcheck2))
//			{}
//			else if(al.contains(forcheck2+"_"+forcheck1))
//			{}
//			else
//			{
//				jishu++;
//				
//				String instName1="inst#gene#unknown#"+forcheck1.substring(5);
//				if(net.nodeExist(instName1))
//				{
//					
//				}
//				else
//				{
//					if(net.nodeExist(forcheck1))
//					{
//						Node node1=new Node(instName1,num1,num2);
//						node1.setType("gene");
//						node1.setPathway("unknown");
//						Node abst1=net.getByName(forcheck1);
//						net.addNode(node1);
//						node1.setBelong(abst1);
//						abst1.addContain(node1);
//						node1.addDown(abst1);
//						node1.addUP(abst1);
//						abst1.addDown(node1);
//						abst1.addUP(node1);
//						Interaction interab=new Interaction(node1,abst1);
//						interab.setType("type");
//						net.addInteraction(interab);
//						interab.setWeight(0.0);
//						
//						interab=new Interaction(abst1, node1);
//						interab.setType("type");
//						net.addInteraction(interab);
//						interab.setWeight(0.0);
//					}
//					else
//					{
//						Node node1=new Node(instName1,num1,num2);
//						node1.setType("gene");
//						node1.setPathway("unknown");
//						Node abst1=new Node(forcheck1,num1,num2);
//						abst1.setType("abst");
//						net.addNode(abst1);
//						net.addNode(node1);
//						node1.setBelong(abst1);
//						abst1.addContain(node1);
//						node1.addDown(abst1);
//						node1.addUP(abst1);
//						abst1.addDown(node1);
//						abst1.addUP(node1);
//						Interaction interab=new Interaction(node1,abst1);
//						interab.setType("type");
//						net.addInteraction(interab);
//						interab.setWeight(0.0);
//						
//						interab=new Interaction(abst1, node1);
//						interab.setType("type");
//						net.addInteraction(interab);
//						interab.setWeight(0.0);
//					}
//				}
//				String instName2="inst#gene#unknown#"+forcheck2.substring(5);
//				if(net.nodeExist(instName2))
//				{
//					
//				}
//				else
//				{
//					if(net.nodeExist(forcheck2))
//					{
//						Node node2=new Node(instName2,num1,num2);
//						node2.setType("gene");
//						node2.setPathway("unknown");
//						Node abst2=net.getByName(forcheck2);
//						net.addNode(node2);
//						node2.setBelong(abst2);
//						abst2.addContain(node2);
//						node2.addDown(abst2);
//						node2.addUP(abst2);
//						abst2.addDown(node2);
//						abst2.addUP(node2);
//						Interaction interab=new Interaction(node2,abst2);
//						interab.setType("type");
//						net.addInteraction(interab);
//						interab.setWeight(0.0);
//						
//						interab=new Interaction(abst2, node2);
//						interab.setType("type");
//						net.addInteraction(interab);
//						interab.setWeight(0.0);
//					}
//					else
//					{
//						Node node2=new Node(instName2,num1,num2);
//						node2.setType("gene");
//						node2.setPathway("unknown");
//						Node abst2=new Node(forcheck2,num1,num2);
//						abst2.setType("abst");
//						net.addNode(abst2);
//						net.addNode(node2);
//						node2.setBelong(abst2);
//						abst2.addContain(node2);
//						node2.addDown(abst2);
//						node2.addUP(abst2);
//						abst2.addDown(node2);
//						abst2.addUP(node2);
//						Interaction interab=new Interaction(node2,abst2);
//						interab.setType("type");
//						net.addInteraction(interab);
//						interab.setWeight(0.0);
//						
//						interab=new Interaction(abst2, node2);
//						interab.setType("type");
//						net.addInteraction(interab);
//						interab.setWeight(0.0);
//					}
//				}
//				Node nd1=net.getByName(instName1);
//				Node nd2=net.getByName(instName2);
//				Interaction interab=new Interaction(nd1,nd2);
//				interab.setType("PPI");
//				net.addInteraction(interab);
//				interab.setWeight(9.0);
//				
//				interab=new Interaction(nd2, nd1);
//				interab.setType("PPI");
//				net.addInteraction(interab);
//				interab.setWeight(9.0);
//				
//				nd1.addDown(nd2);
//				nd1.addUP(nd2);
//				nd2.addDown(nd1);
//				nd2.addUP(nd1);
//				
//				al.add(forcheck1+"_"+forcheck2);
//				
//			}
//			
//		}
//		System.out.println("# of ppi edge after DIP: "+jishu);
		
		Hashtable<String, String> ht_StringmapABST=new Hashtable<String, String>();
		//	String forcheck="";
			
		File Stringmapfile=new File(Stringmap);
		BufferedReader reader = new BufferedReader(new FileReader(Stringmapfile));
		String line="";
		while ((line = reader.readLine()) != null)
		{
			String[] str=line.split("\t");
			if(str.length<2)
				continue;
			
			else
			{
				String StringID=str[0];
				String keggid="abst#"+str[1].split(";")[0];   //seperated by ";"
				ht_StringmapABST.put(StringID, keggid);
			}
			
		}
		
		File StringNET=new File(Stringname);
		reader = new BufferedReader(new FileReader(StringNET));
//		reader.readLine();
		while ((line = reader.readLine()) != null)
		{
//			System.out.println("IO: "+line);
			String[] str=line.split(" ");
//			System.out.println(line);
//			if(Integer.valueOf(str[2])<900)
//				continue;
			String forcheck1="";
			String forcheck2="";
			
				String id1=str[0];
				if(ht_StringmapABST.containsKey(id1))
				{
					forcheck1=ht_StringmapABST.get(id1);
				}
				else
				{
					forcheck1="abst#"+id1;
				}
				
			
			
				String id2=str[1];
				if(ht_StringmapABST.containsKey(id2))
				{
					forcheck2=ht_StringmapABST.get(id2);
				}
				else
				{
					forcheck2="abst#"+id2;
				}
			
			if(al.contains(forcheck1+"_"+forcheck2))
			{}
			else if(al.contains(forcheck2+"_"+forcheck1))
			{}
			else
			{
				jishu++;
				String instName1="inst#gene#unknown#"+forcheck1.substring(5);
				if(net.nodeExist(instName1))
				{
					
				}
				else
				{
					if(net.nodeExist(forcheck1))
					{
						Node node1=new Node(instName1,num1,num2);
						node1.setType("gene");
						node1.setPathway("unknown");
						Node abst1=net.getByName(forcheck1);
						net.addNode(node1);
						node1.setBelong(abst1);
						abst1.addContain(node1);
						node1.addDown(abst1);
						node1.addUP(abst1);
						abst1.addDown(node1);
						abst1.addUP(node1);
						Interaction interab=new Interaction(node1,abst1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst1, node1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
					else
					{
						Node node1=new Node(instName1,num1,num2);
						node1.setType("gene");
						node1.setPathway("unknown");
						Node abst1=new Node(forcheck1,num1,num2);
						abst1.setType("abst");
						net.addNode(abst1);
						net.addNode(node1);
						node1.setBelong(abst1);
						abst1.addContain(node1);
						node1.addDown(abst1);
						node1.addUP(abst1);
						abst1.addDown(node1);
						abst1.addUP(node1);
						Interaction interab=new Interaction(node1,abst1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst1, node1);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
				}
				String instName2="inst#gene#unknown#"+forcheck2.substring(5);
				if(net.nodeExist(instName2))
				{
					
				}
				else
				{
					if(net.nodeExist(forcheck2))
					{
						Node node2=new Node(instName2,num1,num2);
						node2.setType("gene");
						node2.setPathway("unknown");
						Node abst2=net.getByName(forcheck2);
						net.addNode(node2);
						node2.setBelong(abst2);
						abst2.addContain(node2);
						node2.addDown(abst2);
						node2.addUP(abst2);
						abst2.addDown(node2);
						abst2.addUP(node2);
						Interaction interab=new Interaction(node2,abst2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst2, node2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
					else
					{
						Node node2=new Node(instName2,num1,num2);
						node2.setType("gene");
						node2.setPathway("unknown");
						Node abst2=new Node(forcheck2,num1,num2);
						abst2.setType("abst");
						net.addNode(abst2);
						net.addNode(node2);
						node2.setBelong(abst2);
						abst2.addContain(node2);
						node2.addDown(abst2);
						node2.addUP(abst2);
						abst2.addDown(node2);
						abst2.addUP(node2);
						Interaction interab=new Interaction(node2,abst2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
						
						interab=new Interaction(abst2, node2);
						interab.setType("type");
						net.addInteraction(interab);
						interab.setWeight(0.0);
					}
				}
				Node nd1=net.getByName(instName1);
				Node nd2=net.getByName(instName2);
				Interaction interab=new Interaction(nd1,nd2);
				interab.setType("PPI");
				net.addInteraction(interab);
				interab.setWeight(9.0);
				
				interab=new Interaction(nd2, nd1);
				interab.setType("PPI");
				net.addInteraction(interab);
				interab.setWeight(9.0);
				
				nd1.addDown(nd2);
				nd1.addUP(nd2);
				nd2.addDown(nd1);
				nd2.addUP(nd1);
				
				al.add(forcheck1+"_"+forcheck2);
				
			}
			
		}
		
//		while(true)
//		{
//			
//		}
		System.out.println("# of PPI after String: "+jishu);
		System.out.println("# of node after String: "+net.nodes.size());
		return net;
	}
	
	
	
	public Network readNetworkfromFile(String filename, int num1, int num2) throws IOException
	{
		int jishu=0; ///////////////count
		
		Network network=new Network();
		
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		int index=0;
		String tempString = null;
    	while ((tempString = reader.readLine()) != null)
    	{
    		
    		String[] care=tempString.split("\t");			
			Node node1=new Node(care[0],num1,num2);

			Node node2=new Node(care[2],num1,num2);

			String edgetype=care[1];
			
			if(network.nodeExist(node1.getname()))
			{
				node1=network.getByName(node1.getname());
			}
			else
			{
				String pattern = "inst#([^#]*)#([^#]+)#[^#]*";
        		Pattern r=Pattern.compile(pattern);
        		Matcher m=r.matcher(node1.getname());
        		node1.setType("abst");
        		if(m.find())
        		{  			
        			node1.setPathway(m.group(2));
        			node1.setType(m.group(1));
        		}       		
				node1.setIndex(index++);
				network.addNode(node1);
			}
			
			if(network.nodeExist(node2.getname()))
			{
				node2=network.getByName(node2.getname());
			}
			else
			{
				String pattern = "inst#([^#]*)#([^#]+)#[^#]*";
        		Pattern r=Pattern.compile(pattern);
        		Matcher m=r.matcher(node2.getname());
        		node2.setType("abst");
        		if(m.find())
        		{  			
        			node2.setPathway(m.group(2));
        			node2.setType(m.group(1));
        		}       		
				node2.setIndex(index++);
				network.addNode(node2);
			}
			
			Interaction interab=new Interaction(node1,node2);
			if(network.edgeExist(node1.getname()+node2.getname()))
			{
				String endnames=node1.getname()+node2.getname();
				interab=network.getByendnames(endnames);
				if(!interab.getType().equals(edgetype))
				{
					String newtype=interab.getType()+"#"+edgetype;
					interab.setType(newtype);
				}
			}
			else
			{
				if(edgetype.equals("type"))
				{
					interab.setType(edgetype);
					network.addInteraction(interab);
					interab.setWeight(0.0);
					
					interab=new Interaction(node2, node1);
					interab.setType(edgetype);
					network.addInteraction(interab);
					interab.setWeight(0.0);
					
					node1.addDown(node2);
					node1.addUP(node2);
					node2.addDown(node1);
					node2.addUP(node1);
					
					node1.setBelong(node2);
					node2.addContain(node1);
				}
				if(edgetype.equals("compound"))
				{
					interab.setType(edgetype);
					network.addInteraction(interab);
					interab.setWeight(1.0);
					
					//unmark from here
//					Node node3=new Node(care[3],num1,num2);
//					if(network.nodeExist(node3.getname()))
//					{
//						node3=network.getByName(node3.getname());
//					}
//					else
//					{
//						node3.setType("compound");
//						network.addNode(node3);
//					}
//					interab.setCompound(node3);
					//unmark ends here
					
					jishu++;
					
					node1.addDown(node2);
					node2.addUP(node1);
				}
				else 
				{
					interab.setType(edgetype);
					network.addInteraction(interab);
					interab.setWeight(1.0);
					
					jishu++;
					
					node1.addDown(node2);
					node2.addUP(node1);
				}
			}
			
    	}
    	reader.close();
    	System.out.println("# of edge: "+jishu);
    	System.out.println("# of node: "+network.nodes.size());
    	//omit the matrix construction
    	return network;
	}
	
	public Network readNetworkfordream(String filename, int num1, int num2) throws IOException
	{
		
		Network network=new Network();
		
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		int index=0;
		String tempString = null;
    	while ((tempString = reader.readLine()) != null)
    	{
    		
    		String[] care=tempString.split("\t");			
			Node node1=new Node(care[0],num1,num2);

			Node node2=new Node(care[2],num1,num2);

			String edgetype=care[1];
			double edgevalue=Double.valueOf(care[3]);
			if(network.nodeExist(node1.getname()))
			{
				node1=network.getByName(node1.getname());
			}
			else
			{
				String pattern = "inst#([^#]*)#([^#]+)#[^#]*";
        		Pattern r=Pattern.compile(pattern);
        		Matcher m=r.matcher(node1.getname());
        		node1.setType("abst");
        		if(m.find())
        		{  			
        			node1.setPathway(m.group(2));
        			node1.setType(m.group(1));
        		}       		
				node1.setIndex(index++);
				network.addNode(node1);
			}
			
			if(network.nodeExist(node2.getname()))
			{
				node2=network.getByName(node2.getname());
			}
			else
			{
				String pattern = "inst#([^#]*)#([^#]+)#[^#]*";
        		Pattern r=Pattern.compile(pattern);
        		Matcher m=r.matcher(node2.getname());
        		node2.setType("abst");
        		if(m.find())
        		{  			
        			node2.setPathway(m.group(2));
        			node2.setType(m.group(1));
        		}       		
				node2.setIndex(index++);
				network.addNode(node2);
			}
			
			Interaction interab=new Interaction(node1,node2);
			if(network.edgeExist(node1.getname()+node2.getname()))
			{
				String endnames=node1.getname()+node2.getname();
				interab=network.getByendnames(endnames);
				if(!interab.getType().equals(edgetype))
				{
					String newtype=interab.getType()+"#"+edgetype;
					interab.setType(newtype);
				}
			}
			else
			{
				if(edgetype.equals("type"))
				{
					interab.setType(edgetype);
					network.addInteraction(interab);
					interab.setWeight(0.0);
					
					interab=new Interaction(node2, node1);
					interab.setType(edgetype);
					network.addInteraction(interab);
					interab.setWeight(0.0);
					
					node1.addDown(node2);
					node1.addUP(node2);
					node2.addDown(node1);
					node2.addUP(node1);
					
					node1.setBelong(node2);
					node2.addContain(node1);
				}
				if(edgetype.equals("compound"))
				{
					interab.setType(edgetype);
					network.addInteraction(interab);
					interab.setWeight(1.0);
					
//					Node node3=new Node(care[3],num1,num2);
//					if(network.nodeExist(node3.getname()))
//					{
//						node3=network.getByName(node3.getname());
//					}
//					else
//					{
//						node3.setType("compound");
//						network.addNode(node3);
//					}
//					interab.setCompound(node3);
					
					node1.addDown(node2);
					node2.addUP(node1);
				}
				else 
				{
					interab.setType(edgetype);
					network.addInteraction(interab);
					interab.setWeight(edgevalue);
					
					node1.addDown(node2);
					node2.addUP(node1);
				}
			}
			
    	}
    	reader.close();
    	//omit the matrix construction
    	return network;
	}
	
	public ArrayList<Node> readConfidVectorfromFile(String cfilename,Network net) throws IOException
	{
		File cfile=new File(cfilename);
		BufferedReader creader = new BufferedReader(new FileReader(cfile));
		ArrayList<Node> al = new ArrayList<Node>();
		
    	String tempString = null;
		while ((tempString = creader.readLine()) != null) 
		{
			String[] tem=tempString.split("\t");
			Node conf=net.getByName(tem[0]);
			if(conf!=null)
			{
				al.add(conf);
				if(tem.length==1)
				{
					conf.setWeightflag(10.0);
					for(Node inst : conf.getContain())
		    		{
		 //   			inst.setFlag(inst.getFlag()+"C");
		    			inst.setWeightflag(10.0);                    //here is a parameter
		    		}
				}
				if(tem.length==2)
				{
					conf.setWeightflag(Double.valueOf(tem[1]));
					for(Node inst : conf.getContain())
		    		{
		 //   			inst.setFlag(inst.getFlag()+"C");
		    			inst.setWeightflag(Double.valueOf(tem[1]));            //here is a parameter
		    		}
				}
	    		
			}
    		
		}
		creader.close();
		return al;
	}
	
	public ArrayList<Node> readStartVectorfromFile(String filename,Network net) throws IOException
	{
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		ArrayList<Node> al = new ArrayList<Node>();
		
		String tempString;
    	while ((tempString = reader.readLine()) != null) 
    	{	
    		if(net.nodeExist(tempString))
    		{
    			Node start=net.getByName(tempString);
 
        		al.add(start);
        		for(Node inst : start.getContain())
        		{
     //   			inst.setFlag(inst.getFlag()+"S");
        		}
    		}
    		else
    		{
    			System.out.println("Start gene "+tempString+" does not exist");
    		}
    			
    	}
    	reader.close();
    	return al;
	}
	
	public ArrayList<Node> readendVectorfromFile(String filename,Network net) throws IOException
	{
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		ArrayList<Node> al = new ArrayList<Node>();
		
		String tempString;
    	while ((tempString = reader.readLine()) != null) 
    	{	
    		if(net.nodeExist(tempString))
    		{
    			Node start=net.getByName(tempString);
    			start.setIsdif(true);
        		al.add(start);
        		for(Node inst : start.getContain())
        		{
        			start.setIsdif(true);
  //      			inst.setFlag(inst.getFlag()+"E");
        		}
    		}
    		else
    		{
    			System.out.println("End gene "+tempString+" does not exist");
    		}
    			
    	}
    	reader.close();
    	return al;
	}
	public ArrayList<Node> readNendVectorfromFile(String filename,Network net,int N) throws IOException
	{
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		ArrayList<Node> al = new ArrayList<Node>();
		int k=0;
		
		String tempString;
    	while ((tempString = reader.readLine()) != null) 
    	{	
    		if(net.nodeExist(tempString))
    		{
    			k++;
    			Node start=net.getByName(tempString);
    			start.setIsdif(true);
        		al.add(start);
        		for(Node inst : start.getContain())
        		{
        			start.setIsdif(true);
  //      			inst.setFlag(inst.getFlag()+"E");
        		}
    		}
    		else
    		{
    			System.out.println("End gene "+tempString+" does not exist");
    		}
    		if(k==N)
    			break;
    			
    	}
    	reader.close();
    	return al;
	}
	
	
	
	public ArrayList<Node> edgeWeight_CuffDiff(String filename, Network net, double percent) throws IOException
	{
//		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
//		{
//			Node nd=i.getValue();
//			if(nd.getType().equals("abst"))
//			{
//				nd.setWeight(0.001);
//
//			}
//		}
		
		ArrayList<Node> target=new ArrayList<Node>();
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		while ((tempString = reader.readLine()) != null) 
		{
			
			String[] str=tempString.split("\t");
			String ndname=str[2];
			if(net.nodeExist(ndname))
			{
		//		System.out.println("IO: "+tempString);
				if((!str[9].equals("inf"))&&(!str[9].equals("#NAME?"))&&(!str[9].equals("-inf")))
				{
			//		System.out.println("1");
					
					Node nd=net.getByName(ndname);
					double d=Double.valueOf(str[9]);
					nd.setFc(d);
					double w=1/(Math.abs(d)+0.001)/nd.getWeightflag();
					nd.setWeight(w);
					
					if((Math.abs(d)>1)&&(str[13].equals("yes")))
					{
						
						target.add(nd);
					}
				}
			}	
		}
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if(nd.getType().equals("abst"))
			{
					for(Node inst : nd.getContain())
		    		{
						inst.setFc(nd.getFc());
						inst.setWeight(nd.getWeight());
		    		}

			}
		}
		
		System.out.println("IO.java: "+target.size());
		
		ArrayList<Node> result=this.sortIncrement(target);
		ArrayList<Node> newtarget=new ArrayList<Node>();
		int count=(int) Math.floor(result.size()*percent);
		
		
		for(int i=0;i<count;i++)
		{
			newtarget.add(result.get(i));
//			System.out.println(result.get(i).getname()+" "+result.get(i).getFc());
		}
		
		System.out.println(newtarget.size());
		return newtarget;
		
	}
	public void edgeWeight_dfgiven2(String filename, Network net) throws IOException
	{
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		int exp_num=0;
		while ((tempString = reader.readLine()) != null) 
		{
			
			String[] str=tempString.split("\t");
			String ndname=str[0];
			
			if(net.nodeExist(ndname))
			{
				
				Node nd=net.getByName(ndname);
				exp_num++;
				double d=Double.valueOf(str[1]);
				nd.setFc(d);
				d=Double.valueOf(str[2]);
				if(d!=1.0)
					nd.setPvalue(d);
				double w=1/-trans(nd.getPvalue(),2);
				nd.setWeight(w/nd.getWeightflag());
					
				for(Node inst : nd.getContain())
	    		{
					inst.setFc(nd.getFc());
					inst.setPvalue(nd.getPvalue());
					inst.setWeight(nd.getWeight());
	    		}
				 
			}	
		}
		System.out.println("gene with exp_num: "+exp_num);
		net.setExp_num(exp_num);
		
		
		reader.close();

	}
	
	public ArrayList<Node> edgeWeight_dfgiven(String filename, Network net,int target, Pathway P) throws IOException
	{
		ArrayList<Node> endlist=new ArrayList<Node>();
		ArrayList<Node> alllist=new ArrayList<Node>();
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		int exp_num=0;
		while ((tempString = reader.readLine()) != null) 
		{
			
			String[] str=tempString.split("\t");
			String ndname=str[0];
			
			if(net.nodeExist(ndname))
			{
				
				Node nd=net.getByName(ndname);
				exp_num++;
				alllist.add(nd);
				double d=Double.valueOf(str[1]);
				nd.setFc(d);
				d=Double.valueOf(str[2]);
				if(d!=1.0)
					nd.setPvalue(d);
				double w=1/-trans(nd.getPvalue(),2);
				nd.setWeight(w/nd.getWeightflag());
				if(nd.getPvalue()<=0.01)
				{
					endlist.add(nd);
					alllist.add(nd);
				}
				else
				{
					alllist.add(nd);
				}
					
				for(Node inst : nd.getContain())
	    		{
					inst.setFc(nd.getFc());
					inst.setPvalue(nd.getPvalue());
					inst.setWeight(nd.getWeight());
	    		}
				 
			}	
		}
		System.out.println("gene with exp_num: "+exp_num);
		net.setExp_num(exp_num);
		
		
		
		P.getABCD(net, endlist);
		P.node_list=endlist;
		
		reader.close();
		System.out.println("IO: target size is: "+endlist.size());
		
		if(endlist.size()>target)
		{
			ArrayList<Node> result=this.sortIncrement(endlist);
			ArrayList<Node> newtarget=new ArrayList<Node>();
	//		int count=(int) Math.floor(result.size()*percent);
			
			
			for(int i=0;i<target;i++)
			{
				newtarget.add(result.get(i));
//				System.out.println(result.get(i).getname()+" "+result.get(i).getFc());
			}
			
			System.out.println("IO: new target size is: "+newtarget.size());
			return newtarget;
		}
		if(endlist.size()<target)
		{
			ArrayList<Node> result=this.sortIncrement(alllist);
			ArrayList<Node> newtarget=new ArrayList<Node>();
		
			for(int i=0;i<target;i++)
			{
				newtarget.add(result.get(i));
		//		System.out.println("=========================="+result.get(i).getWeight());
			}
					
			System.out.println("IO: new target size is: "+newtarget.size());
			return newtarget;
		}
		
		return endlist;
	}
	
	public ArrayList<Node> edgeWeight_noexp(Network net, int target, Pathway P) throws IOException
	{
		ArrayList<Node> endlist=new ArrayList<Node>();
		ArrayList<Node> alllist=new ArrayList<Node>();
		net.setExp_num(0);
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if(nd.getType().equals("abst"))
			{

				if(Math.abs(nd.getFc())>1)
				{
					endlist.add(nd);
					alllist.add(nd);
					for(Node inst : nd.getContain())
		    		{
						inst.setFc(nd.getFc());
						inst.setWeight(nd.getWeight());
		    		}
				}
				else
				{
					alllist.add(nd);
					for(Node inst : nd.getContain())
		    		{
						inst.setFc(nd.getFc());
						inst.setWeight(nd.getWeight());
		    		}
				}

//					for(Node inst : nd.getContain())
//		    		{
//						inst.setFc(nd.getFc());
//						inst.setWeight(nd.getWeight());
//		    		}
					

	

			}
		}
		P.getABCD(net, endlist);
		P.node_list=endlist;
		System.out.println("IO: target size is: "+endlist.size());
		
		if(endlist.size()>target)
		{
			ArrayList<Node> result=this.sortIncrement(endlist);
			ArrayList<Node> newtarget=new ArrayList<Node>();
	//		int count=(int) Math.floor(result.size()*percent);
			
			
			for(int i=0;i<target;i++)
			{
				newtarget.add(result.get(i));
//				System.out.println(result.get(i).getname()+" "+result.get(i).getFc());
			}
			
			System.out.println("IO: new target size is: "+newtarget.size());
			return newtarget;
		}
		if(endlist.size()<target)
		{
			ArrayList<Node> result=this.sortIncrement(alllist);
			ArrayList<Node> newtarget=new ArrayList<Node>();
		
			for(int i=0;i<target;i++)
			{
				newtarget.add(result.get(i));
		//		System.out.println("=========================="+result.get(i).getWeight());
			}
					
			System.out.println("IO: new target size is: "+newtarget.size());
			return newtarget;
		}
		return endlist;
		
	}
	
	public ArrayList<Node> edgeWeight_NOreplicate(String filename, Network net, int target, Pathway P) throws IOException
	{
		ArrayList<Node> endlist=new ArrayList<Node>();
		ArrayList<Node> alllist=new ArrayList<Node>();
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		
		int exp_num=0;
		while ((tempString = reader.readLine()) != null) 
		{
			String[] str=tempString.split("\t");
			String ndname=str[0];
			
			if(net.nodeExist(ndname))
			{
				Node nd=net.getByName(ndname);
				
//				nd.setUnisample(Double.valueOf(str[1]));
//				nd.setUnisample(Double.valueOf(str[2]));
				double d1=Double.valueOf(str[1]);
				double d2=Double.valueOf(str[2]);
				double[] v1={d1};
				double[] v2={d2};
				nd.setExpressionfromV(v1, 1);
				nd.setExpressionfromV(v2, 2);
				
				exp_num++;
				nd.setFc(trans(d2/d1,2));
				if(Math.abs(d2/d1-1)<=0.01)
				{
					
					nd.setWeight(69);//1/log2(1.01)
				}
				else
				{
					
					double w=1/Math.abs(nd.getFc());
					nd.setWeight(w/nd.getWeightflag());
				}
				
		        
			}	
		}
		System.out.println("gene with exp_num: "+exp_num);
		net.setExp_num(exp_num);
		
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if(nd.getType().equals("abst"))
			{

				
					if(Math.abs(nd.getFc())>1)
					{
						endlist.add(nd);
						alllist.add(nd);
						for(Node inst : nd.getContain())
			    		{
		//					inst.setFlag(inst.getFlag()+"E");
			    			
							
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    		}
					}
					else
					{
						alllist.add(nd);
						for(Node inst : nd.getContain())
			    		{
		//					inst.setFlag(inst.getFlag()+"E");
			    			
							
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    		}
					}
				

					
					
					////////////////////////////////////////////////////////////////////////////////	
//					double v=-trans(nd.getPvalue(),2);
//				other.write(nd.getname()+"\tGene\t"+v);	
//				other.write("\r\n");
					
				/////////////////////////////////////////////////////////////////////////////////////
	

			}
		}
		
		P.getABCD(net, endlist);
		P.node_list=endlist;
		
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if(nd.getType().equals("group"))
			{
				nd.setExpression4group(1, 1, net);
				nd.setExpression4group(2, 1, net);
				double d1=nd.getExpression(1)[0];
				double d2=nd.getExpression(2)[0];
			
				
				nd.setFc(trans(d2/d1,2));
				if(Math.abs(d2/d1-1)<=0.01)
				{
					
					nd.setWeight(69);//1/log2(1.01)
				}
				else
				{
					
					double w=1/Math.abs(nd.getFc());
					nd.setWeight(w/nd.getWeightflag());
				}
				
			
			
		
		        
			}
		}
		reader.close();
		System.out.println("IO: target size is: "+endlist.size());
		
		if(endlist.size()>target)
		{
			ArrayList<Node> result=this.sortIncrement(endlist);
			ArrayList<Node> newtarget=new ArrayList<Node>();
	//		int count=(int) Math.floor(result.size()*percent);
			
			
			for(int i=0;i<target;i++)
			{
				newtarget.add(result.get(i));
//				System.out.println(result.get(i).getname()+" "+result.get(i).getFc());
			}
			
			System.out.println("IO: new target size is: "+newtarget.size());
			return newtarget;
		}
		if(endlist.size()<target)
		{
			ArrayList<Node> result=this.sortIncrement(alllist);
			ArrayList<Node> newtarget=new ArrayList<Node>();
		
			for(int i=0;i<target;i++)
			{
				newtarget.add(result.get(i));
		//		System.out.println("=========================="+result.get(i).getWeight());
			}
					
			System.out.println("IO: new target size is: "+newtarget.size());
			return newtarget;
		}
		
		
		return endlist;
	}
	
	public ArrayList<Node> edgeWeight4Inte(String filename_mrna,String filename_pro, Network net, int num1, int num2, int target, Pathway P) throws IOException
	{
		ArrayList<Node> endlist_mrna=new ArrayList<Node>();
		ArrayList<Node> endlist_pro=new ArrayList<Node>();
		ArrayList<Node> alllist_mrna=new ArrayList<Node>();
		ArrayList<Node> alllist_pro=new ArrayList<Node>();
		File file=new File(filename_mrna);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		
			int exp_num_mrna=0;
			while ((tempString = reader.readLine()) != null) 
			{
				
				String[] str=tempString.split("\t");
				if(str.length!=1+num1+num2)
					continue;
				String ndname=str[0];
	//			System.out.println(ndname+" "+str.length);
				if(net.nodeExist(ndname))
				{
					
					Node nd=net.getByName(ndname);
					
					
					String[] con1=Arrays.copyOfRange(str, 1, 1+num1);
					String[] con2=Arrays.copyOfRange(str, 1+num1, 1+num1+num2);
					nd.setExpressionfromarray(con1, 1);
					nd.setExpressionfromarray(con2, 2);
					
					
					exp_num_mrna++;
					alllist_mrna.add(nd);
					//cal node fold change
					DescriptiveStatistics stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(1).length; i++) {
			            stats.addValue(nd.getExpression(1)[i]);
			            }
			        double mean1=stats.getMean();
			
			        stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(2).length; i++) {
			            stats.addValue(nd.getExpression(2)[i]);
			            }
			        double mean2=stats.getMean();
			
	//		        System.out.println(nd.getname()+" "+mean1+" "+mean2);
			        if(!((mean1==0)||(mean2==0)))
			        {
			        	double dev=mean2/mean1;
				        double fc=trans(dev, 2);
				        nd.setFc(fc);
			        }
			        
					
			        double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
			        nd.setPvalue(s);
			        
			
				}	
			}
			System.out.println("mrna with exp_num: "+exp_num_mrna);
			
			file=new File(filename_pro);
			reader = new BufferedReader(new FileReader(file));
			
				int exp_num_pro=0;
				while ((tempString = reader.readLine()) != null) 
				{
					
					String[] str=tempString.split("\t");
					if(str.length!=1+num1+num2)
						continue;
					String ndname=str[0];
		//			System.out.println(ndname+" "+str.length);
					if(net.nodeExist(ndname))
					{
						
						Node nd=net.getByName(ndname);
						
						
						String[] con1=Arrays.copyOfRange(str, 1, 1+num1);
						String[] con2=Arrays.copyOfRange(str, 1+num1, 1+num1+num2);
						nd.setExpressionfromarray(con1, 1);
						nd.setExpressionfromarray(con2, 2);
						
						
						exp_num_pro++;
						alllist_pro.add(nd);
						//cal node fold change
						DescriptiveStatistics stats=new DescriptiveStatistics();
				        for( int i = 0; i < nd.getExpression(1).length; i++) {
				            stats.addValue(nd.getExpression(1)[i]);
				            }
				        double mean1=stats.getMean();
				
				        stats=new DescriptiveStatistics();
				        for( int i = 0; i < nd.getExpression(2).length; i++) {
				            stats.addValue(nd.getExpression(2)[i]);
				            }
				        double mean2=stats.getMean();
				
		//		        System.out.println(nd.getname()+" "+mean1+" "+mean2);
				        if(!((mean1==0)||(mean2==0)))
				        {
				        	double dev=mean2/mean1;
					        double fc=trans(dev, 2);
					        nd.setFc(fc);
				        }
				        
						
				        double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
				        nd.setPvalue(s);
				        
				
					}	
				}
				System.out.println("pro with exp_num: "+exp_num_pro);
			
				net.setExp_num(exp_num_mrna+exp_num_pro);
			/////////////////////////////////////////////////////////////////////////////////////
//			FileWriter other = new FileWriter("data/yeast/cell_wall/forNetres.txt");
			///////////////////////////////////////////////////////////////////////////////////////
			for(Map.Entry<String, Node> i : net.getNodes().entrySet())
			{
				
				
				Node nd=i.getValue();
				if(nd.getType().equals("abst"))
				{
					
					double w=1/-trans(nd.getPvalue(),2);
					nd.setWeight(w/nd.getWeightflag());
	
					boolean tt=new TTest().tTest(nd.getExpression(1), nd.getExpression(2), 0.01/net.getExp_num());
	//				boolean tt=new TTest().tTest(nd.getExpression(1), nd.getExpression(2), 0.5);
	//				if(tt && nd.getFc()>=1)
					
					if(tt )	
					{
						if(nd.getname().split("\\.").length>1)
						{
							endlist_pro.add(nd);
			//				alllist_pro.add(nd);
						}
						else
						{
							endlist_mrna.add(nd);
			//				alllist_mrna.add(nd);
						}
						
						
						for(Node inst : nd.getContain())
			    		{
		//					inst.setFlag(inst.getFlag()+"E");
			    			
							inst.setPvalue(nd.getPvalue());
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    		}
						
						
						////////////////////////////////////////////////////////////////////////////////	
	//					double v=-trans(nd.getPvalue(),2);
	//				other.write(nd.getname()+"\tGene\t"+v);	
	//				other.write("\r\n");
						
					/////////////////////////////////////////////////////////////////////////////////////
					}
					else
					{
						
						for(Node inst : nd.getContain())
			    		{
							inst.setPvalue(nd.getPvalue());
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    		}
						/////////////////////////////////////////////////////////////////////////////////////		
	//					double v=-trans(nd.getPvalue(),2);
	//				other.write(nd.getname()+"\tGene\t"+v);	
	//				other.write("\r\n");
					/////////////////////////////////////////////////////////////////////////////////////
					}

				}
			}
			ArrayList<Node> endlist=new ArrayList<Node>();
			endlist.addAll(endlist_pro);
			endlist.addAll(endlist_mrna);
			
			
			P.getABCD(net, endlist);
			P.node_list=endlist;
			///////////
	//		other.close();
			/////////////
			for(Map.Entry<String, Node> i : net.getNodes().entrySet())
			{
				Node nd=i.getValue();
				if(nd.getType().equals("group"))
				{
					nd.setExpression4group(1, num1, net);
					nd.setExpression4group(2, num2, net);
					double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
					nd.setPvalue(s);
					double w=1/-trans(nd.getPvalue(),2);
					nd.setWeight(w/nd.getWeightflag());
					
//					DescriptiveStatistics stats=new DescriptiveStatistics();
//			        for( int k = 0; k < nd.getExpression(1).length; k++) {
//			            stats.addValue(nd.getExpression(1)[k]);
//			            }
//			        double mean1=stats.getMean();
//			        stats=new DescriptiveStatistics();
//			        for( int k = 0; k < nd.getExpression(2).length; k++) {
//			            stats.addValue(nd.getExpression(2)[k]);
//			            }
//			        double mean2=stats.getMean();
//			        double fc=Math.abs(mean1/mean2);
//			        nd.setFc(fc);
			        
				}
			}
		
		
//			for(Map.Entry<String, Interaction> i : net.getInteractions().entrySet())
//			{
//				Interaction inter= i.getValue();
//				if(inter.getType().equals("type"))
//				{
//					inter.setWeight(0);
//				}
//				else
//				{
//					Node nd1=inter.getNode1();
//					Node nd2=inter.getNode2();
//					double w=(nd1.getFc()+nd2.getFc())/2;
//					inter.setWeight(w);
//				}
//			}
			
			
			System.out.println("IO: target size is: "+endlist.size());
			System.out.println("IO: target mrna size is: "+endlist_mrna.size());
			System.out.println("IO: target pro size is: "+endlist_pro.size());
			System.out.println("IO: alllist pro size is: "+alllist_pro.size());
			System.out.println("IO: alllist mrna size is: "+alllist_mrna.size());
			reader.close();
			
			ArrayList<Node> newtarget=new ArrayList<Node>();
			if(endlist_pro.size()>target)
			{
				ArrayList<Node> result=this.sortIncrement(endlist_pro);
				for(int i=0;i<target;i++)
				{
					newtarget.add(result.get(i));
				}
				System.out.println("IO: new target size after pro is: "+newtarget.size());
				
			}
			else
			{
				ArrayList<Node> result=this.sortIncrement(alllist_pro);
				for(int i=0;i<target;i++)
				{
					newtarget.add(result.get(i));
				}
				System.out.println("IO: new target size after pro is: "+newtarget.size());
				
			}
			
			
	//		newtarget.addAll(endlist_pro);
			if(endlist_mrna.size()>target)
			{
				ArrayList<Node> result=this.sortIncrement(endlist_mrna);
				for(int i=0;i<target;i++)
				{
					newtarget.add(result.get(i));
				}
				System.out.println("IO: new target size is: "+newtarget.size());
				return newtarget;
			}
			if(endlist_mrna.size()<target)
			{
				ArrayList<Node> result=this.sortIncrement(alllist_mrna);
				for(int i=0;i<target;i++)
				{
					newtarget.add(result.get(i));
				}
				System.out.println("IO: new target size is: "+newtarget.size());
				return newtarget;
			}
				
			
		
		return endlist;
		
	}
	
	public ArrayList<Node> edgeWeight4metab(String filename_mrna, String filename_metab, Network net, int num1, int num2, int target, Pathway P,int num1meta, int num2meta) throws IOException
	{
		File file=new File(filename_metab);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		
			int exp_num=0;
			while ((tempString = reader.readLine()) != null) 
			{
				
				String[] str=tempString.split("\t");
				if(str.length!=1+num1meta+num2meta)
					continue;
				String ndname=str[0];
				if(net.nodeExist(ndname))
				{
					
					Node nd=net.getByName(ndname);
					
					
					String[] con1=Arrays.copyOfRange(str, 1, 1+num1meta);
					String[] con2=Arrays.copyOfRange(str, 1+num1meta, 1+num1meta+num2meta);
					nd.setExpressionfromarray(con1, 1);
					nd.setExpressionfromarray(con2, 2);
					
					
					//cal node fold change
					DescriptiveStatistics stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(1).length; i++) {
			            stats.addValue(nd.getExpression(1)[i]);
			            }
			        double mean1=stats.getMean();
			
			        stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(2).length; i++) {
			            stats.addValue(nd.getExpression(2)[i]);
			            }
			        double mean2=stats.getMean();
			
	//		        System.out.println(nd.getname()+" "+mean1+" "+mean2);
			        if(!((mean1==0)||(mean2==0)))
			        {
			        	double dev=mean2/mean1;
				        double fc=trans(dev, 2);
				        nd.setFc(fc);
				        double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
				        if(s!=1.0)
				        	nd.setPvalue(s);
			        }
					
					double w=1/-trans(nd.getPvalue(),2);
					nd.setWeight(w/nd.getWeightflag());
					
					for(Node inst : nd.getContain())
			    	{	
							inst.setPvalue(nd.getPvalue());
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    	}
					  
				}	
			}
		reader.close();
		
		for(Map.Entry<String, Interaction> i : net.getInteractions().entrySet())
		{
			Interaction inter=i.getValue();
			Node nd=inter.getCompound();
			inter.setWeight(nd.getWeight());
		}
		
		
		
		ArrayList<Node> endlist=new ArrayList<Node>();
		ArrayList<Node> alllist=new ArrayList<Node>();
		file=new File(filename_mrna);
		reader = new BufferedReader(new FileReader(file));
		
			exp_num=0;
			while ((tempString = reader.readLine()) != null) 
			{
				
				String[] str=tempString.split("\t");
				if(str.length!=1+num1+num2)
					continue;
				String ndname=str[0];
	//			System.out.println(ndname+" "+str.length);
				if(net.nodeExist(ndname))
				{
					
					Node nd=net.getByName(ndname);
					
					
					String[] con1=Arrays.copyOfRange(str, 1, 1+num1);
					String[] con2=Arrays.copyOfRange(str, 1+num1, 1+num1+num2);
					nd.setExpressionfromarray(con1, 1);
					nd.setExpressionfromarray(con2, 2);
					
					
					exp_num++;
					//cal node fold change
					DescriptiveStatistics stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(1).length; i++) {
			            stats.addValue(nd.getExpression(1)[i]);
			            }
			        double mean1=stats.getMean();
			
			        stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(2).length; i++) {
			            stats.addValue(nd.getExpression(2)[i]);
			            }
			        double mean2=stats.getMean();
			
	//		        System.out.println(nd.getname()+" "+mean1+" "+mean2);
			        if(!((mean1==0)||(mean2==0)))
			        {
			        	double dev=mean2/mean1;
				        double fc=trans(dev, 2);
				        nd.setFc(fc);
				        double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
				        if(s!=1.0)
				        	nd.setPvalue(s);
			        }
			        
					
			        
			        
			
				}	
			}
			System.out.println("gene with exp_num: "+exp_num);
			net.setExp_num(exp_num);
			/////////////////////////////////////////////////////////////////////////////////////
//			FileWriter other = new FileWriter("data/yeast/cell_wall/forNetres.txt");
			///////////////////////////////////////////////////////////////////////////////////////
			for(Map.Entry<String, Node> i : net.getNodes().entrySet())
			{
				
				
				Node nd=i.getValue();
				String pattern = "cpd:";
        		Pattern r=Pattern.compile(pattern);
        		Matcher m=r.matcher(nd.getname());
        		if(m.find())
        			continue;
				if(nd.getType().equals("abst"))
				{
					
					double w=1/-trans(nd.getPvalue(),2);
					nd.setWeight(w/nd.getWeightflag());
	//				System.out.println(nd.getname()+"\t"+nd.getPvalue()+"\t"+nd.getWeight());
					boolean tt=new TTest().tTest(nd.getExpression(1), nd.getExpression(2), 0.01/net.getExp_num());
	//				boolean tt=new TTest().tTest(nd.getExpression(1), nd.getExpression(2), 0.5);
	//				if(tt && nd.getFc()>=1)
					
					if(tt )	
					{
						endlist.add(nd);
						alllist.add(nd);
						for(Node inst : nd.getContain())
			    		{
		//					inst.setFlag(inst.getFlag()+"E");
			    			
							inst.setPvalue(nd.getPvalue());
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    		}
						
						
						////////////////////////////////////////////////////////////////////////////////	
	//					double v=-trans(nd.getPvalue(),2);
	//				other.write(nd.getname()+"\tGene\t"+v);	
	//				other.write("\r\n");
						
					/////////////////////////////////////////////////////////////////////////////////////
					}
					else
					{
						alllist.add(nd);
						for(Node inst : nd.getContain())
			    		{
							inst.setPvalue(nd.getPvalue());
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    		}
						/////////////////////////////////////////////////////////////////////////////////////		
	//					double v=-trans(nd.getPvalue(),2);
	//				other.write(nd.getname()+"\tGene\t"+v);	
	//				other.write("\r\n");
					/////////////////////////////////////////////////////////////////////////////////////
					}

				}
			}
			P.getABCD(net, endlist);
			P.node_list=endlist;
			///////////
	//		other.close();
			/////////////
			for(Map.Entry<String, Node> i : net.getNodes().entrySet())
			{
				Node nd=i.getValue();
				if(nd.getType().equals("group"))
				{
					nd.setExpression4group(1, num1, net);
					nd.setExpression4group(2, num2, net);
					double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
					nd.setPvalue(s);
					double w=1/-trans(nd.getPvalue(),2);
					nd.setWeight(w/nd.getWeightflag());
					
//					DescriptiveStatistics stats=new DescriptiveStatistics();
//			        for( int k = 0; k < nd.getExpression(1).length; k++) {
//			            stats.addValue(nd.getExpression(1)[k]);
//			            }
//			        double mean1=stats.getMean();
//			        stats=new DescriptiveStatistics();
//			        for( int k = 0; k < nd.getExpression(2).length; k++) {
//			            stats.addValue(nd.getExpression(2)[k]);
//			            }
//			        double mean2=stats.getMean();
//			        double fc=Math.abs(mean1/mean2);
//			        nd.setFc(fc);
			        
				}
			}
		
		
//			for(Map.Entry<String, Interaction> i : net.getInteractions().entrySet())
//			{
//				Interaction inter= i.getValue();
//				if(inter.getType().equals("type"))
//				{
//					inter.setWeight(0);
//				}
//				else
//				{
//					Node nd1=inter.getNode1();
//					Node nd2=inter.getNode2();
//					double w=(nd1.getFc()+nd2.getFc())/2;
//					inter.setWeight(w);
//				}
//			}
			System.out.println("IO: target size is: "+endlist.size());
			if(endlist.size()>target)
			{
				ArrayList<Node> result=this.sortIncrement(endlist);
				ArrayList<Node> newtarget=new ArrayList<Node>();
		//		int count=(int) Math.floor(result.size()*percent);
				
				
				for(int i=0;i<target;i++)
				{
					newtarget.add(result.get(i));
//					System.out.println(result.get(i).getname()+" "+result.get(i).getFc());
				}
				
				System.out.println("IO: new target size is: "+newtarget.size());
				return newtarget;
			}
			if(endlist.size()<target)
			{
				ArrayList<Node> result=this.sortIncrement(alllist);
				ArrayList<Node> newtarget=new ArrayList<Node>();
			
				for(int i=0;i<target;i++)
				{
					newtarget.add(result.get(i));
			//		System.out.println("=========================="+result.get(i).getWeight());
				}
						
				System.out.println("IO: new target size is: "+newtarget.size());
				return newtarget;
			}
			
			

		reader.close();
		return endlist;
		
	}
	
	public ArrayList<Node> edgeWeight(String filename, Network net, int num1, int num2, int target, Pathway P) throws IOException
	{
		ArrayList<Node> endlist=new ArrayList<Node>();
		ArrayList<Node> alllist=new ArrayList<Node>();
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		
			int exp_num=0;
			while ((tempString = reader.readLine()) != null) 
			{
				
				String[] str=tempString.split("\t");
				if(str.length!=1+num1+num2)
					continue;
				String ndname=str[0];
	//			System.out.println(ndname+" "+str.length);
				if(net.nodeExist(ndname))
				{
					
					Node nd=net.getByName(ndname);
					
					
					String[] con1=Arrays.copyOfRange(str, 1, 1+num1);
					String[] con2=Arrays.copyOfRange(str, 1+num1, 1+num1+num2);
					nd.setExpressionfromarray(con1, 1);
					nd.setExpressionfromarray(con2, 2);
					
					
					exp_num++;
					//cal node fold change
					DescriptiveStatistics stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(1).length; i++) {
			            stats.addValue(nd.getExpression(1)[i]);
			            }
			        double mean1=stats.getMean();
			
			        stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(2).length; i++) {
			            stats.addValue(nd.getExpression(2)[i]);
			            }
			        double mean2=stats.getMean();
			
	//		        System.out.println(nd.getname()+" "+mean1+" "+mean2);
			        if(!((mean1==0)||(mean2==0)))
			        {
			        	double dev=mean2/mean1;
				        double fc=trans(dev, 2);
				        nd.setFc(fc);
				        double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
				        if(s!=1.0)
				        	nd.setPvalue(s);
			        }
			        
					
			        
			        
			
				}	
			}
			System.out.println("gene with exp_num: "+exp_num);
			net.setExp_num(exp_num);
			/////////////////////////////////////////////////////////////////////////////////////
//			FileWriter other = new FileWriter("data/yeast/cell_wall/forNetres.txt");
			///////////////////////////////////////////////////////////////////////////////////////
			for(Map.Entry<String, Node> i : net.getNodes().entrySet())
			{
				
				
				Node nd=i.getValue();
				if(nd.getType().equals("abst"))
				{
					
					double w=1/-trans(nd.getPvalue(),2);
					nd.setWeight(w/nd.getWeightflag());
	//				System.out.println(nd.getname()+"\t"+nd.getPvalue()+"\t"+nd.getWeight());
					boolean tt=new TTest().tTest(nd.getExpression(1), nd.getExpression(2), 0.01/net.getExp_num());
	//				boolean tt=new TTest().tTest(nd.getExpression(1), nd.getExpression(2), 0.5);
	//				if(tt && nd.getFc()>=1)
					
					if(tt )	
					{
						endlist.add(nd);
						alllist.add(nd);
						for(Node inst : nd.getContain())
			    		{
		//					inst.setFlag(inst.getFlag()+"E");
			    			
							inst.setPvalue(nd.getPvalue());
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    		}
						
						
						////////////////////////////////////////////////////////////////////////////////	
	//					double v=-trans(nd.getPvalue(),2);
	//				other.write(nd.getname()+"\tGene\t"+v);	
	//				other.write("\r\n");
						
					/////////////////////////////////////////////////////////////////////////////////////
					}
					else
					{
						alllist.add(nd);
						for(Node inst : nd.getContain())
			    		{
							inst.setPvalue(nd.getPvalue());
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    		}
						/////////////////////////////////////////////////////////////////////////////////////		
	//					double v=-trans(nd.getPvalue(),2);
	//				other.write(nd.getname()+"\tGene\t"+v);	
	//				other.write("\r\n");
					/////////////////////////////////////////////////////////////////////////////////////
					}

				}
			}
			P.getABCD(net, endlist);
			P.node_list=endlist;
			///////////
	//		other.close();
			/////////////
			for(Map.Entry<String, Node> i : net.getNodes().entrySet())
			{
				Node nd=i.getValue();
				if(nd.getType().equals("group"))
				{
					nd.setExpression4group(1, num1, net);
					nd.setExpression4group(2, num2, net);
					double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
					nd.setPvalue(s);
					double w=1/-trans(nd.getPvalue(),2);
					nd.setWeight(w/nd.getWeightflag());
					
//					DescriptiveStatistics stats=new DescriptiveStatistics();
//			        for( int k = 0; k < nd.getExpression(1).length; k++) {
//			            stats.addValue(nd.getExpression(1)[k]);
//			            }
//			        double mean1=stats.getMean();
//			        stats=new DescriptiveStatistics();
//			        for( int k = 0; k < nd.getExpression(2).length; k++) {
//			            stats.addValue(nd.getExpression(2)[k]);
//			            }
//			        double mean2=stats.getMean();
//			        double fc=Math.abs(mean1/mean2);
//			        nd.setFc(fc);
			        
				}
			}
		
		
//			for(Map.Entry<String, Interaction> i : net.getInteractions().entrySet())
//			{
//				Interaction inter= i.getValue();
//				if(inter.getType().equals("type"))
//				{
//					inter.setWeight(0);
//				}
//				else
//				{
//					Node nd1=inter.getNode1();
//					Node nd2=inter.getNode2();
//					double w=(nd1.getFc()+nd2.getFc())/2;
//					inter.setWeight(w);
//				}
//			}
			System.out.println("IO: target size is: "+endlist.size());
			if(endlist.size()>target)
			{
				ArrayList<Node> result=this.sortIncrement(endlist);
				ArrayList<Node> newtarget=new ArrayList<Node>();
		//		int count=(int) Math.floor(result.size()*percent);
				
				
				for(int i=0;i<target;i++)
				{
					newtarget.add(result.get(i));
//					System.out.println(result.get(i).getname()+" "+result.get(i).getFc());
				}
				
				System.out.println("IO: new target size is: "+newtarget.size());
				return newtarget;
			}
			if(endlist.size()<target)
			{
				ArrayList<Node> result=this.sortIncrement(alllist);
				ArrayList<Node> newtarget=new ArrayList<Node>();
			
				for(int i=0;i<target;i++)
				{
					newtarget.add(result.get(i));
			//		System.out.println("=========================="+result.get(i).getWeight());
				}
						
				System.out.println("IO: new target size is: "+newtarget.size());
				return newtarget;
			}
			
			

		reader.close();
		return endlist;
		
	}
	public void edgeWeight2(String filename, Network net, int num1, int num2) throws IOException
	{
		
		File file=new File(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		
			int exp_num=0;
			while ((tempString = reader.readLine()) != null) 
			{
				
				String[] str=tempString.split("\t");
				if(str.length!=1+num1+num2)
					continue;
				String ndname=str[0];
	//			System.out.println(ndname+" "+str.length);
				if(net.nodeExist(ndname))
				{
					
					Node nd=net.getByName(ndname);
					
					
					String[] con1=Arrays.copyOfRange(str, 1, 1+num1);
					String[] con2=Arrays.copyOfRange(str, 1+num1, 1+num1+num2);
					nd.setExpressionfromarray(con1, 1);
					nd.setExpressionfromarray(con2, 2);
					
					
					exp_num++;
					//cal node fold change
					DescriptiveStatistics stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(1).length; i++) {
			            stats.addValue(nd.getExpression(1)[i]);
			            }
			        double mean1=stats.getMean();
			
			        stats=new DescriptiveStatistics();
			        for( int i = 0; i < nd.getExpression(2).length; i++) {
			            stats.addValue(nd.getExpression(2)[i]);
			            }
			        double mean2=stats.getMean();
			
	//		        System.out.println(nd.getname()+" "+mean1+" "+mean2);
			        if(!((mean1==0)||(mean2==0)))
			        {
			        	double dev=mean2/mean1;
				        double fc=trans(dev, 2);
				        nd.setFc(fc);
			        }
			        
					
			        double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
			        nd.setPvalue(s);
			        
			
				}	
			}
			System.out.println("gene with exp_num: "+exp_num);
			net.setExp_num(exp_num);
			/////////////////////////////////////////////////////////////////////////////////////
//			FileWriter other = new FileWriter("data/yeast/cell_wall/forNetres.txt");
			///////////////////////////////////////////////////////////////////////////////////////
			for(Map.Entry<String, Node> i : net.getNodes().entrySet())
			{
				
				
				Node nd=i.getValue();
				if(nd.getType().equals("abst"))
				{
					
					double w=1/-trans(nd.getPvalue(),2);
					nd.setWeight(w/nd.getWeightflag());
	
					boolean tt=new TTest().tTest(nd.getExpression(1), nd.getExpression(2), 0.01/net.getExp_num());
	//				boolean tt=new TTest().tTest(nd.getExpression(1), nd.getExpression(2), 0.5);
	//				if(tt && nd.getFc()>=1)
	
						for(Node inst : nd.getContain())
			    		{
							inst.setPvalue(nd.getPvalue());
			    			inst.setExpressionfromV(nd.getExpression(1), 1);
							inst.setExpressionfromV(nd.getExpression(2), 2);
							inst.setFc(nd.getFc());
		//					inst.setWeight(inst.getFc()*inst.getWeightflag());
							inst.setWeight(nd.getWeight());
			    		}
						/////////////////////////////////////////////////////////////////////////////////////		
	//					double v=-trans(nd.getPvalue(),2);
	//				other.write(nd.getname()+"\tGene\t"+v);	
	//				other.write("\r\n");
					/////////////////////////////////////////////////////////////////////////////////////
					

				}
			}
			///////////
	//		other.close();
			/////////////
			for(Map.Entry<String, Node> i : net.getNodes().entrySet())
			{
				Node nd=i.getValue();
				if(nd.getType().equals("group"))
				{
					nd.setExpression4group(1, num1, net);
					nd.setExpression4group(2, num2, net);
					double s=new TTest().tTest(nd.getExpression(1), nd.getExpression(2));
					nd.setPvalue(s);
					double w=1/-trans(nd.getPvalue(),2);
					nd.setWeight(w/nd.getWeightflag());
					
//					DescriptiveStatistics stats=new DescriptiveStatistics();
//			        for( int k = 0; k < nd.getExpression(1).length; k++) {
//			            stats.addValue(nd.getExpression(1)[k]);
//			            }
//			        double mean1=stats.getMean();
//			        stats=new DescriptiveStatistics();
//			        for( int k = 0; k < nd.getExpression(2).length; k++) {
//			            stats.addValue(nd.getExpression(2)[k]);
//			            }
//			        double mean2=stats.getMean();
//			        double fc=Math.abs(mean1/mean2);
//			        nd.setFc(fc);
			        
				}
			}
		
		
//			for(Map.Entry<String, Interaction> i : net.getInteractions().entrySet())
//			{
//				Interaction inter= i.getValue();
//				if(inter.getType().equals("type"))
//				{
//					inter.setWeight(0);
//				}
//				else
//				{
//					Node nd1=inter.getNode1();
//					Node nd2=inter.getNode2();
//					double w=(nd1.getFc()+nd2.getFc())/2;
//					inter.setWeight(w);
//				}
//			}

			

		reader.close();
		
		
	}
	
	/////////////////////////////////////////////dream_time///////////////////////////////////////////////////
	public ArrayList<Node> edgeWeight_dream(String dir, Network net, int num1) throws IOException
	{
		ArrayList<Node> endlist=new ArrayList<Node>();
		
		File file=new File(dir);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		
		while ((tempString = reader.readLine()) != null) 
		{
			String[] str=tempString.split("\t");
			if(str.length!=1+num1)
				continue;
			String ndname=str[0];
			
			if(net.nodeExist(ndname))
			{
				Node nd=net.getByName(ndname);	
				endlist.add(nd);
				String[] con1=Arrays.copyOfRange(str, 1, 1+num1);
				nd.setExpressionfromarray(con1, 1);
				for(Node inst : nd.getContain())
		    	{
		    		inst.setExpressionfromV(nd.getExpression(1), 1);
		    	}
			}	
		}
		
	
		
		
		
		return endlist;
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	public void edgeWeight_time(String dir, Network net, int num1) throws IOException
	{
		File file=new File(dir);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		
		while ((tempString = reader.readLine()) != null) 
		{
			String[] str=tempString.split("\t");
			if(str.length!=1+num1)
				continue;
			String ndname=str[0];
			
			if(net.nodeExist(ndname))
			{
				Node nd=net.getByName(ndname);		
				String[] con1=Arrays.copyOfRange(str, 1, 1+num1);
				nd.setExpressionfromarray(con1, 1);
			}	
		}
		
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if(nd.getType().equals("abst"))
			{
				
				for(Node inst : nd.getContain())
		    	{
		    		inst.setExpressionfromV(nd.getExpression(1), 1);
		    	}
			}
	
		}
		
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if(nd.getType().equals("group"))
			{
				nd.setExpression4group(1, num1, net);
			}
		}
		
		for(Map.Entry<String, Interaction> i : net.getInteractions().entrySet())
		{
			Interaction inter= i.getValue();
			if(inter.getType().equals("type"))
			{
				inter.setWeight(0);
			}
			else
			{
				Node nd1=inter.getNode1();
				Node nd2=inter.getNode2();
				double[] n1=nd1.getExpression(1);
				double[] n2=nd2.getExpression(1);
				Variance v=new Variance();
				if((v.evaluate(n1)==0)||(v.evaluate(n2)==0))
				{
					inter.setWeight(1*inter.getWeight());
				}
				else
				{
					double k=new PearsonsCorrelation().correlation(n1,n2);
					k=1-Math.abs(k);
					inter.setWeight(k*inter.getWeight());
				}

				
			}
		}
		reader.close();
	}
	public String edgeWeight_time2(String expdir, Network net, int tnum, int repnum) throws IOException, REngineException, REXPMismatchException
	{
		this.edgeWeight_time(expdir, net, tnum*repnum);
		RConnection c = new RConnection();
		REXP x=new REXP();
		int[] k={tnum,repnum};
		String a=new File("").getAbsoluteFile().getAbsolutePath();
		System.out.println(a);
//		c.assign("expdir", "C:/Users/install/Documents/EclipseWorkplace/KEGG_multi_dynamic/"+expdir);
		c.assign("expdir", a+"/"+expdir);
		String dir=a+"/"+expdir;
		int occ=dir.lastIndexOf("/");
		String relwd=dir.substring(0, occ);
		System.out.println(relwd);
		
		c.assign("paras", k);
//		c.assign("relwd", "C:/Users/install/Documents/EclipseWorkplace/KEGG_multi_dynamic/"+relwd);
		c.assign("relwd", relwd);
		
		
		File file=new File("dif_maSigPro_dev.R");
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String tempString;
		while((tempString=reader.readLine())!=null)
		{
			System.out.println(tempString);
			REXP xx=c.eval(tempString);
			if(tempString.length()>=3)
			{
				if(tempString.substring(0, 3).equals("sum"))
				{
					System.out.println("satisfy threshold num: "+xx.asInteger());
				}
			}
			
			
			
		}
	//	String[] xx=x.asStrings();
	//	System.out.println(xx[0]);
	//	System.out.println(xx[1]);
		
		c.close();
		return relwd;
		
		
	}
	
	public ArrayList<Node> sortIncrement(ArrayList<Node> al)
	{
		Collections.sort(al,new Comparator<Node>()
		{
			public int compare(Node arg0, Node arg1) 
			{
				Double darg0=arg0.getWeight();
				Double darg1=arg1.getWeight();
	//			double tem=(arg0.getWeight() - arg1.getWeight())*100000;
	//			return (int)tem;
				return darg0.compareTo(darg1);
            }
        });
		return al;
	}
	
}
