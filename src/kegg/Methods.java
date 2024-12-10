package kegg;
import java.io.*;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class Methods 
{
	//store instances information, function: updateArraylist, writeArrayList
	public ArrayList<Node> al=new ArrayList<Node>();
	//store abst information, function: updateHash, writeHast, 2 different situations: for PN or for MN
	//these two are for PN
	public Hashtable<String, Integer> htOutdegree=new Hashtable<String, Integer>();
	public Hashtable<String, Integer> htIndegree=new Hashtable<String, Integer>();
	//these two are for MN
	public Hashtable<String, Integer> htOutdegree_MN=new Hashtable<String, Integer>();
	public Hashtable<String, Integer> htIndegree_MN=new Hashtable<String, Integer>();
	//this is for PN
	public String relation;
	//this is for MN
	public String reaction;
		
	//output abst in Hashtable, these abst have 0 outdegree or indegree
	public void writeHash(Hashtable<String, Integer> ht, FileWriter fw) throws IOException
	{
		for(String str : ht.keySet())
		{
			if(ht.get(str)==0)
			{
				fw.write(str);
				fw.write("\r\n");
			}
		}
	}	
	//record abst indegree and outdegree
	public void updateHash(Hashtable<String, Integer> out, Hashtable<String, Integer> in, String abst, Character cha)
	{
		if(cha.equals('+'))
		{
			if(out.containsKey(abst))
			{
				int newnum=out.get(abst)+1;
				out.put(abst, newnum);
			}
			else 
			{
				out.put(abst, 1);
			}
			
			if(in.containsKey(abst))
			{
				//do nothing
			}
			else 
			{
				in.put(abst, 0);
			}
		}
		if(cha.equals('-'))
		{
			if(in.containsKey(abst))
			{
				int newnum=in.get(abst)+1;
				in.put(abst, newnum);
			}
			else
			{
				in.put(abst, 1);
			}
			
			if(out.containsKey(abst))
			{
				//do nothing
			}
			else 
			{
				out.put(abst, 0);
			}
		}
	}
	//output nodes in the arraylist, these nodes have 0 outdegree or indegree
	public void writeArraylist(ArrayList<Node> al,FileWriter fw)
	{
	
	}
	//record nodes indegree and outdegree, cha("+" or "-") indicates it is a source or end
	public void updateArraylist(ArrayList<Node> al,String fullname,Character cha)
	{
		if(cha.equals('+'))
		{
			int counter=0;
			for(Node a : al)
			{
				counter++;
				if(a.getName().equals(fullname))
				{
					a.setOutdegree(a.getOutdegree()+1);
					break;
				}
				if(counter==al.size()) 
				{
					Node newAdd=new Node(fullname);
					newAdd.setOutdegree(1);
					al.add(newAdd);
				}	
			}
		}
		if(cha.equals('-'))
		{
			int counter=0;
			for(Node a : al)
			{
				if(a.getName().equals(fullname))
				{
					a.setIndegree(a.getIndegree()+1);
					break;
				}
				if(counter==al.size()) 
				{
					Node newAdd=new Node(fullname);
					newAdd.setIndegree(1);
					al.add(newAdd);
				}
			}
			
		}	
	}
	//find pathways files
	public File[] getPathfiles(String dir)
	{
		File file=new File(dir);
		File[] tempList = file.listFiles();
		return tempList;
	}
	
	//entry is a gene
	public void metGene(Hashtable<String,String> ht1,Hashtable<String,String> ht2,String id,String name,String type)
	{
		ht1.put(id, name.toString().replaceAll(" ", "_"));
    	ht2.put(id, type);
	}
	
	//entry is a compound
	public void metCompound(Hashtable<String,String> ht1,Hashtable<String,String> ht2,String id,String name,String type)
	{
		ht1.put(id, name.toString().replaceAll(" ", "_"));
    	ht2.put(id, type);
	}
	
	//entry is a map
	public void metMap(Hashtable<String,String> ht1,Hashtable<String,String> ht2,String id,String name,String type)
	{
		ht1.put(id, name.toString().replaceAll(" ", "_"));
    	ht2.put(id, type);
	}
	
	//entry is a group
	public void metGroup(Hashtable<String,String> ht1,Hashtable<String,String> ht2,String id,String name,String type,BufferedReader reader) throws IOException
	{
		String nextString=null;
		Pattern pattern=Pattern.compile("</entry>");
		String groupname="|com|";
		while(!(pattern.matcher(nextString=reader.readLine()).find()))
    	{
			Pattern pattern11=Pattern.compile("<component id=\"(.*)\"");
	        Matcher m11=pattern11.matcher(nextString);
	        if(m11.find())
	        {
	        	groupname=groupname+ht1.get(m11.group(1))+"|com|";
	        }
    	}
		int count=groupname.length()-5;
		groupname=groupname.substring(0, count);
		ht1.put(id, groupname.replaceAll(" ", "_"));
		ht2.put(id, type);
	}
	
	//record group components information
	public void recGroup(StringBuilder sb, String groupname, Hashtable<String, Integer> ht3, String pname)
	{
		String[] str_out=groupname.split("\\|com\\|");
		for(int j=1;j<str_out.length;j++)
		{
			String[] str_in=str_out[j].split("_");
			{
				for(int i=0;i<str_in.length;i++)
				{
					String inst="inst#gene#"+pname+"#"+str_in[i];
					
					if(!ht3.containsKey(inst))
					{
						sb.append(inst+"\t"+"type"+"\t"+"abst#"+str_in[i]+"\r\n");
						ht3.put(inst, 1);
					}
				}
			}
		}
	}
	
	//reaction function
	public String react(Boolean ifrever, Hashtable<String,Integer> ht3, BufferedReader reader,String pname) throws IOException
	{
		String type="";
		if(ifrever)
		{
			type="reversible";
		}
		else {
			type="irreversible";
		}
		StringBuilder reactString = new StringBuilder();
		ArrayList<String> substrate=new ArrayList<String>();
		ArrayList<String> product=new ArrayList<String>();
		
		String nextString=null;
		Pattern pattern=Pattern.compile("</reaction>");
		while(!(pattern.matcher(nextString=reader.readLine()).find()))
		{
			Pattern pattern31=Pattern.compile("<substrate id=\"(.*)\" name=\"(.*)\"");
			Matcher m31=pattern31.matcher(nextString);
			Pattern pattern32=Pattern.compile("<product id=\"(.*)\" name=\"(.*)\"");
			Matcher m32=pattern32.matcher(nextString);
			if(m31.find())
			{
				String cid=m31.group(1);
				String rname=m31.group(2);
				String cname="inst#compound#"+pname+"#"+m31.group(2);
				if(!ht3.containsKey(cname))
				{
					reactString.append(cname+"\t"+"type"+"\t"+"abst#"+rname+"\r\n");
					ht3.put(cname, 1);
				}
				substrate.add(rname);
			}
			if(m32.find())
			{
				String cid=m32.group(1);
				String rname=m32.group(2);
				String cname="inst#compound#"+pname+"#"+m32.group(2);
				if(!ht3.containsKey(cname))
				{
					reactString.append(cname+"\t"+"type"+"\t"+"abst#"+rname+"\r\n");
					ht3.put(cname, 1);
				}
				product.add(rname);
			}
		}
		if(substrate.size()==0||product.size()==0)
		{
			return "";
		}
		else {
			for(int i=0;i<substrate.size();i++)
			{
				for(int j=0;j<product.size();j++)
				{
					reactString.append("inst#compound#"+pname+"#"+substrate.get(i)+"\t"+type+"\t"+"inst#compound#"+pname+"#"+product.get(j)+"\r\n");
					updateHash(this.htOutdegree_MN, this.htIndegree_MN, "abst#"+substrate.get(i), '+');
					updateHash(this.htOutdegree_MN, this.htIndegree_MN, "abst#"+product.get(j), '-');
					if(type.equals("reversible"))
					{
						reactString.append("inst#compound#"+pname+"#"+product.get(j)+"\t"+type+"\t"+"inst#compound#"+pname+"#"+substrate.get(i)+"\r\n");
						updateHash(this.htOutdegree_MN, this.htIndegree_MN, "abst#"+substrate.get(i), '-');
						updateHash(this.htOutdegree_MN, this.htIndegree_MN, "abst#"+product.get(j), '+');
					}
				}
			}
		}
		return reactString.toString();
	}
	
	//relation is nongroup--nongroup
	public String relatNN(String id1, String id2,Hashtable<String,String> ht1, Hashtable<String,String> ht2, Hashtable<String,Integer> ht3, BufferedReader reader,String pname) throws IOException
	{
		StringBuilder relationString = new StringBuilder();
		String[] set1=ht1.get(id1).split("_");
		String[] set2=ht1.get(id2).split("_");
		String type1=ht2.get(id1);
		String type2=ht2.get(id2);
		
		String nextString=null;
		Pattern pattern=Pattern.compile("</relation>");
		while(!(pattern.matcher(nextString=reader.readLine()).find()))
    	{
			Pattern pattern21=Pattern.compile("<subtype name=\"(.*)\" value=\"(.*)\"");
	        Matcher m21=pattern21.matcher(nextString);
	        if(m21.find())
	        {
	        	String subtype=m21.group(1);
	        	
	        	String cid=m21.group(2);
	        	String cname="";
	        	Boolean index=false;
	        	if(subtype.equals("compound"))
				{
					index=true;
					cname="inst#compound#"+pname+"#"+ht1.get(cid);
					if(!ht3.containsKey(cname))
					{
						relationString.append(cname+"\t"+"type"+"\t"+"abst#"+ht1.get(cid)+"\r\n");
    					ht3.put(cname, 1);
					}
				}
	        	
	        	for(int i=0;i<set1.length;i++)
	    		{
	        		String fullname1="inst#"+type1+"#"+pname+"#"+set1[i];
	        		if(!ht3.containsKey(fullname1))
    				{
    					relationString.append(fullname1+"\t"+"type"+"\t"+"abst#"+set1[i]+"\r\n");
    					ht3.put(fullname1, 1);
    				}
	        		
	       // 		updateArraylist(this.al, fullname1, '+');
	        		
	        		
	    			for(int j=0;j<set2.length;j++)
	    			{	    				
	    				String fullname2="inst#"+type2+"#"+pname+"#"+set2[j];	    				
	    				if(!ht3.containsKey(fullname2))
	    				{
	    					relationString.append(fullname2+"\t"+"type"+"\t"+"abst#"+set2[j]+"\r\n");
	    					ht3.put(fullname2, 1);
	    				}
	    				
	    				if(index)
	    				{
	    					relationString.append(fullname1+"\t"+subtype+"\t"+fullname2+"\t"+cname+"\r\n");
	    				}
	    				else 
	    				{
	    					relationString.append(fullname1+"\t"+subtype+"\t"+fullname2+"\r\n");
						}
	    			    				
	    	//			updateArraylist(this.al, fullname2, '-');	    	
	    				updateHash(this.htOutdegree, this.htIndegree, "abst#"+set1[i], '+');
	    				updateHash(this.htOutdegree, this.htIndegree, "abst#"+set2[j], '-');
	    			}
	    		}

	        }
    	}
		return relationString.toString();
	}
	
	//relation is group--nongroup
	public String relatGN(String id1, String id2,Hashtable<String,String> ht1, Hashtable<String,String> ht2, Hashtable<String,Integer> ht3, BufferedReader reader,String pname) throws IOException
	{
		StringBuilder relationString = new StringBuilder();
		String[] set2=ht1.get(id2).split("_");
		String type1=ht2.get(id1);
		String type2=ht2.get(id2);
		
		String nextString=null;
		Pattern pattern=Pattern.compile("</relation>");
		while(!(pattern.matcher(nextString=reader.readLine()).find()))
    	{
			Pattern pattern21=Pattern.compile("<subtype name=\"(.*)\" value=\"(.*)\"");
	        Matcher m21=pattern21.matcher(nextString);
	        if(m21.find())
	        {
	        	String subtype=m21.group(1);
	        	
	        	String cid=m21.group(2);
	        	String cname="";
	        	Boolean index=false;
	        	if(subtype.equals("compound"))
				{
					index=true;
					cname="inst#compound#"+pname+"#"+ht1.get(cid);
				}
	        	
	        	String fullname1="inst#"+type1+"#"+pname+"#"+ht1.get(id1);
				if(!ht3.containsKey(fullname1))
				{
					relationString.append(fullname1+"\t"+"type"+"\t"+"abst#"+ht1.get(id1)+"\r\n");
					ht3.put(fullname1, 1);
				}
				recGroup(relationString, fullname1, ht3, pname);
		//		updateArraylist(this.al, fullname1, '+');
				
				
	        	for(int j=0;j<set2.length;j++)
    			{
    				String fullname2="inst#"+type2+"#"+pname+"#"+set2[j];
    				if(!ht3.containsKey(fullname2))
    				{
    					relationString.append(fullname2+"\t"+"type"+"\t"+"abst#"+set2[j]+"\r\n");
    					ht3.put(fullname2, 1);
    				}
    				
    				if(index)
    				{
    					relationString.append(fullname1+"\t"+subtype+"\t"+fullname2+"\t"+cname+"\r\n");
    				}
    				else 
    				{
    					relationString.append(fullname1+"\t"+subtype+"\t"+fullname2+"\r\n");
					}
    				
    	//			updateArraylist(this.al, fullname2, '-');
    				updateHash(this.htOutdegree, this.htIndegree, "abst#"+ht1.get(id1), '+');
    				updateHash(this.htOutdegree, this.htIndegree, "abst#"+set2[j], '-');
    			}
	        }
    	}
		return relationString.toString();
	}
	
	//relation is nongroup--group
	public String relatNG(String id1, String id2,Hashtable<String,String> ht1, Hashtable<String,String> ht2, Hashtable<String,Integer> ht3, BufferedReader reader,String pname) throws IOException
	{
		StringBuilder relationString = new StringBuilder();
		String[] set1=ht1.get(id1).split("_");
		String type1=ht2.get(id1);
		String type2=ht2.get(id2);
		
		String nextString=null;
		Pattern pattern=Pattern.compile("</relation>");
		while(!(pattern.matcher(nextString=reader.readLine()).find()))
    	{
			Pattern pattern21=Pattern.compile("<subtype name=\"(.*)\" value=\"(.*)\"");
	        Matcher m21=pattern21.matcher(nextString);
	        if(m21.find())
	        {
	        	String subtype=m21.group(1);
	        	
	        	String cid=m21.group(2);
	        	String cname="";
	        	Boolean index=false;
	        	if(subtype.equals("compound"))
				{
					index=true;
					cname="inst#compound#"+pname+"#"+ht1.get(cid);
				}
	        	
	        	String fullname2="inst#"+type2+"#"+pname+"#"+ht1.get(id2);
				if(!ht3.containsKey(fullname2))
				{
					relationString.append(fullname2+"\t"+"type"+"\t"+"abst#"+ht1.get(id2)+"\r\n");
					ht3.put(fullname2, 1);
				}
				recGroup(relationString, fullname2, ht3, pname);
		//		updateArraylist(this.al, fullname2, '-');
				
				
	        	for(int i=0;i<set1.length;i++)
    			{
    				String fullname1="inst#"+type1+"#"+pname+"#"+set1[i];
    				if(!ht3.containsKey(fullname1))
    				{
    					relationString.append(fullname1+"\t"+"type"+"\t"+"abst#"+set1[i]+"\r\n");
    					ht3.put(fullname1, 1);
    				}
				
    				if(index)
    				{
    					relationString.append(fullname1+"\t"+subtype+"\t"+fullname2+"\t"+cname+"\r\n");
    				}
    				else 
    				{
    					relationString.append(fullname1+"\t"+subtype+"\t"+fullname2+"\r\n");
					}
    				
    	//			updateArraylist(this.al, fullname1, '+');
    				updateHash(this.htOutdegree, this.htIndegree, "abst#"+ht1.get(id2), '-');
    				updateHash(this.htOutdegree, this.htIndegree, "abst#"+set1[i], '+');
    				
    			}
	        }
    	}
		return relationString.toString();
	}
	
	//relation is group--group
	public String relatGG(String id1, String id2,Hashtable<String,String> ht1, Hashtable<String,String> ht2, Hashtable<String,Integer> ht3, BufferedReader reader,String pname) throws IOException
	{
		StringBuilder relationString = new StringBuilder();
		String type1=ht2.get(id1);
		String type2=ht2.get(id2);
		
		String nextString=null;
		Pattern pattern=Pattern.compile("</relation>");
		while(!(pattern.matcher(nextString=reader.readLine()).find()))
    	{
			Pattern pattern21=Pattern.compile("<subtype name=\"(.*)\" value=\"(.*)\"");
	        Matcher m21=pattern21.matcher(nextString);
	        if(m21.find())
	        {
	        	String subtype=m21.group(1);
	        	
	        	String cid=m21.group(2);
	        	String cname="";
	        	Boolean index=false;
	        	if(subtype.equals("compound"))
				{
					index=true;
					cname="inst#compound#"+pname+"#"+ht1.get(cid);
				}
	        	
	        	String fullname1="inst#"+type1+"#"+pname+"#"+ht1.get(id1);
	        	String fullname2="inst#"+type2+"#"+pname+"#"+ht1.get(id2);
	        	if(!ht3.containsKey(fullname1))
				{
					relationString.append(fullname1+"\t"+"type"+"\t"+"abst#"+ht1.get(id1)+"\r\n");
					ht3.put(fullname1, 1);
				}
				if(!ht3.containsKey(fullname2))
				{
					relationString.append(fullname2+"\t"+"type"+"\t"+"abst#"+ht1.get(id2)+"\r\n");
					ht3.put(fullname2, 1);
				}
				recGroup(relationString, fullname1, ht3, pname);
				recGroup(relationString, fullname2, ht3, pname);
				if(index)
				{
					relationString.append(fullname1+"\t"+subtype+"\t"+fullname2+"\t"+cname+"\r\n");
				}
				else 
				{
					relationString.append(fullname1+"\t"+subtype+"\t"+fullname2+"\r\n");
				}
				
		//		updateArraylist(this.al, fullname1, '+');
		//		updateArraylist(this.al, fullname2, '-');				
				updateHash(this.htOutdegree, this.htIndegree, "abst#"+ht1.get(id1), '+');
				updateHash(this.htOutdegree, this.htIndegree, "abst#"+ht1.get(id2), '-');
	        }
    	}
		return relationString.toString();
	}
	
	//read a pathway file, 4PN is for build gene-only network, 4MN is for build metabolism network, 4PMN is for build gene and metabolism network
	public void readPathfile4PPN(String inputdir,String pname) throws IOException
	{
		File rfile=new File(inputdir);
		BufferedReader reader = new BufferedReader(new FileReader(rfile));

		Hashtable<String,String> ht_id2name=new Hashtable<String,String>();    //id-rawname
		Hashtable<String,String> ht_id2type=new Hashtable<String,String>();    //id-type
		Hashtable<String,Integer> ht_abst2rndnum=new Hashtable<String,Integer>();    //inst-1, check if it is already written, for confi gene
		StringBuilder relationString=new StringBuilder();   //expandable variable, used to output results£¬ for PN, coresponding relation
		StringBuilder reactionString=new StringBuilder(); //for MN, coresponding reaction
		
		String tempString = null;
		while ((tempString = reader.readLine()) != null) 
	    {
	    	Pattern pattern1=Pattern.compile("<entry id=\"(.*)\".*name=\"(.*)\" type=\"([^\"]*)\"");
	        Matcher m1=pattern1.matcher(tempString);
	        Pattern pattern2=Pattern.compile("<relation entry1=\"(.*)\" entry2=\"(.*)\" type=\"(.*)\""); 
	        Matcher m2=pattern2.matcher(tempString);
	        Pattern pattern3=Pattern.compile("<reaction id=\"(.*)\" name.*type=\"(.*)\"");
	        Matcher m3=pattern3.matcher(tempString);
	        
	        //match entry
	        if(m1.find())
	        {
	        	String id=m1.group(1);
	        	String name=m1.group(2);
	        	String type=m1.group(3);
	        	if(type.equals("gene"))
	        		metGene(ht_id2name,ht_id2type,id,name,type);
	        	else if(type.equals("compound"))
	        		metCompound(ht_id2name,ht_id2type,id,name,type);
	        	else if(type.equals("map"))
	        	{
	        		//do nothing
	        	}
	        	else if(type.equals("group"))
	        		metGroup(ht_id2name,ht_id2type,id,name,type,reader);
	        }
	        
	        //match relation
	        else if(m2.find())
	        {
	        	String ID1=m2.group(1);
	        	String ID2=m2.group(2);
	        	String type=m2.group(3);
	        	if(type.equals("maplink"))
	        	{
	        		continue;
	        	}
	        	if(!ht_id2name.containsKey(ID1)||!ht_id2name.containsKey(ID2))
	        	{
	        		continue;
	        	}
	        	//in the relation, neither g1 or g2 is group 
	        	if((!ht_id2type.get(ID1).equals("group"))&&(!ht_id2type.get(ID2).equals("group")))
	        	{
	        		
	        		relationString.append(relatNN(ID1, ID2, ht_id2name, ht_id2type, ht_abst2rndnum, reader, pname));  		
	        		
	        	}
	        	//in the relation, g1 is a group, g2 is not
	        	if((ht_id2type.get(ID1).equals("group"))&&(!ht_id2type.get(ID2).equals("group")))
	        	{
	        		relationString.append(relatGN(ID1, ID2, ht_id2name, ht_id2type, ht_abst2rndnum, reader, pname));
	        	}
	        	
	        	//in the relation, g2 is a group, g1 is not
	        	if((ht_id2type.get(ID2).equals("group"))&&(!ht_id2type.get(ID1).equals("group")))
	        	{
	        		relationString.append(relatNG(ID1, ID2, ht_id2name, ht_id2type, ht_abst2rndnum, reader, pname));
	        	}
	        	
	        	//in the relation, both g1 and g2 are group
	        	if((ht_id2type.get(ID1).equals("group"))&&(ht_id2type.get(ID2).equals("group")))
	        	{
	        		relationString.append(relatGG(ID1, ID2, ht_id2name, ht_id2type, ht_abst2rndnum, reader, pname));
	        	}
	        	
	        }
	        
	        //match reaction
	        else if(m3.find())
	        {
	        	String type=m3.group(2);
	        	if(type.equals("irreversible"))
	        	{
	        		reactionString.append(react(false, ht_abst2rndnum, reader, pname));
	        	}
	        	else if(type.equals("reversible"))
	        	{
	        		reactionString.append(react(true, ht_abst2rndnum, reader, pname));
	        	}
	        }
	    }
		//return relationString.toString();
		this.relation=relationString.toString();
		this.reaction=reactionString.toString();
	}
	
	public String getCpdName(BufferedReader reader, Hashtable<String,String> ht1) throws IOException
	{
		String nextString=null;
		Pattern pattern=Pattern.compile("</relation>");
		String cname="";
		while(!(pattern.matcher(nextString=reader.readLine()).find()))
    	{
			Pattern pattern21=Pattern.compile("<subtype name=\"(.*)\" value=\"(.*)\"");
	        Matcher m21=pattern21.matcher(nextString);
	        if(m21.find())
	        {
	        	String subtype=m21.group(1);
	        	
	        	String cid=m21.group(2);
	        	cname="";
	        	if(subtype.equals("compound"))
				{
					cname=ht1.get(cid);
				}
	        }
    	}
		return cname;   	
	}
	public String searchMap(String reqPath, String resFileName,String resPath, String reqGeneName, String cname, String mapPos,Hashtable<String,Integer> ht3) throws IOException
	{
		File rfile=new File(resFileName);
		BufferedReader reader = new BufferedReader(new FileReader(rfile));
		Hashtable<String,String> ht_id2name=new Hashtable<String,String>();    //id-rawname
		Hashtable<String,String> ht_id2type=new Hashtable<String,String>();    //id-type
		StringBuilder relationString=new StringBuilder();   //expandable variable, used to output results£¬ for PN, coresponding relation
		
		String tempString = null;
		while ((tempString = reader.readLine()) != null) 
	    {
	    	Pattern pattern1=Pattern.compile("<entry id=\"(.*)\".*name=\"(.*)\" type=\"([^\"]*)\"");
	        Matcher m1=pattern1.matcher(tempString);
	        Pattern pattern2=Pattern.compile("<relation entry1=\"(.*)\" entry2=\"(.*)\" type=\"(.*)\""); 
	        Matcher m2=pattern2.matcher(tempString);
	        
	        //match entry
	        if(m1.find())
	        {
	        	String id=m1.group(1);
	        	String name=m1.group(2);
	        	String type=m1.group(3);
	        	if(type.equals("gene"))
	        		metGene(ht_id2name,ht_id2type,id,name,type);
	        	else if(type.equals("compound"))
	        		metCompound(ht_id2name,ht_id2type,id,name,type);
	        	else if(type.equals("map"))
	        		metMap(ht_id2name,ht_id2type,id,name,type);
	        }
	        
	        //match relation
	        else if(m2.find())
	        {
	        	String ID1=m2.group(1);
	        	String ID2=m2.group(2);
	        	String type=m2.group(3);
	        	if(!ht_id2name.containsKey(ID1)||!ht_id2name.containsKey(ID2))
	        	{
	        		continue;
	        	}
	        	if(type.equals("maplink"))
	        	{
//	        		System.out.println("checkpoint_Methods.java_2:"+ID1+" "+ID2);
	        		if((ht_id2type.get(ID1).equals("map"))&&(ht_id2type.get(ID2).equals("map")))
	        			continue;
	        		else if((ht_id2type.get(ID1).equals("map"))&&(mapPos.equals("map2")))
//	        		else if((ht_id2type.get(ID1).equals("map")))
	        		{
	        			String temReqPath=ht_id2name.get(ID1);
	        			if (!temReqPath.equals(reqPath))
	        				continue;
	        			String temResGeneName=ht_id2name.get(ID2);
	        			String temCname=this.getCpdName(reader, ht_id2name);
	        			if (!temCname.equals(cname))
	        				continue;
	        			
	        			String[] set1=reqGeneName.split("_");
	        			String[] set2=temResGeneName.split("_");
	        			for(int i=0;i<set1.length;i++)
	    	    		{
	    	        		String fullname1="inst#gene#"+reqPath.split(":")[1]+"#"+set1[i];
	    	        		if(!ht3.containsKey(fullname1))
	        				{
	        					relationString.append(fullname1+"\t"+"type"+"\t"+"abst#"+set1[i]+"\r\n");
	        					ht3.put(fullname1, 1);
	        				}
	    	        		
	    	        		
	    	    			for(int j=0;j<set2.length;j++)
	    	    			{	    				
	    	    				String fullname2="inst#gene#"+resPath.split(":")[1]+"#"+set2[j];	    				
	    	    				if(!ht3.containsKey(fullname2))
	    	    				{
	    	    					relationString.append(fullname2+"\t"+"type"+"\t"+"abst#"+set2[j]+"\r\n");
	    	    					ht3.put(fullname2, 1);
	    	    				}
	    	    				
	    	    				relationString.append(fullname1+"\t"+type+"\t"+fullname2+"\r\n");

	    	    			}
	    	    		}
	        					
	        		}
	        		else if((ht_id2type.get(ID2).equals("map"))&&(mapPos.equals("map1")))
//	        		else if((ht_id2type.get(ID2).equals("map")))
	        		{
	        			String temReqPath=ht_id2name.get(ID2);
	        			if (!temReqPath.equals(reqPath))
	        				continue;
	        			String temResGeneName=ht_id2name.get(ID1);
	        			String temCname=this.getCpdName(reader, ht_id2name);
	        			if (!temCname.equals(cname))
	        				continue;
	        			
	        			String[] set1=temResGeneName.split("_");
	        			String[] set2=reqGeneName.split("_");
	        			for(int i=0;i<set1.length;i++)
	    	    		{
	    	        		String fullname1="inst#gene#"+resPath.split(":")[1]+"#"+set1[i];
	    	        		if(!ht3.containsKey(fullname1))
	        				{
	        					relationString.append(fullname1+"\t"+"type"+"\t"+"abst#"+set1[i]+"\r\n");
	        					ht3.put(fullname1, 1);
	        				}
	    	        		
	    	        		
	    	    			for(int j=0;j<set2.length;j++)
	    	    			{	    				
	    	    				String fullname2="inst#gene#"+reqPath.split(":")[1]+"#"+set2[j];	    				
	    	    				if(!ht3.containsKey(fullname2))
	    	    				{
	    	    					relationString.append(fullname2+"\t"+"type"+"\t"+"abst#"+set2[j]+"\r\n");
	    	    					ht3.put(fullname2, 1);
	    	    				}
	    	    				
	    	    				relationString.append(fullname1+"\t"+type+"\t"+fullname2+"\r\n");

	    	    			}
	    	    		}

	        			
	        		}
	        		
	        	}
	        }
	    }
		return relationString.toString();
	}
	public String readmaplink(String inputfile,String pname,String pathdir) throws IOException
	{
		File rfile=new File(inputfile);
		BufferedReader reader = new BufferedReader(new FileReader(rfile));

		Hashtable<String,String> ht_id2name=new Hashtable<String,String>();    //id-rawname
		Hashtable<String,String> ht_id2type=new Hashtable<String,String>();    //id-type
		Hashtable<String,Integer> ht_abst2rndnum=new Hashtable<String,Integer>();    //inst-1, check if it is already written, for confi gene
		StringBuilder relationString=new StringBuilder();   //expandable variable, used to output results£¬ for PN, coresponding relation
		
		String tempString = null;
		while ((tempString = reader.readLine()) != null) 
	    {
	    	Pattern pattern1=Pattern.compile("<entry id=\"(.*)\".*name=\"(.*)\" type=\"([^\"]*)\"");
	        Matcher m1=pattern1.matcher(tempString);
	        Pattern pattern2=Pattern.compile("<relation entry1=\"(.*)\" entry2=\"(.*)\" type=\"(.*)\""); 
	        Matcher m2=pattern2.matcher(tempString);
	        
	        //match entry
	        if(m1.find())
	        {
	        	String id=m1.group(1);
	        	String name=m1.group(2);
	        	String type=m1.group(3);
	        	if(type.equals("gene"))
	        		metGene(ht_id2name,ht_id2type,id,name,type);
	        	else if(type.equals("compound"))
	        		metCompound(ht_id2name,ht_id2type,id,name,type);
	        	else if(type.equals("map"))
	        		metMap(ht_id2name,ht_id2type,id,name,type);
//	        	System.out.println("checkpoint_Methods.java_1:"+id);
	        }
	        
	        //match relation
	        else if(m2.find())
	        {
	        	String ID1=m2.group(1);
	        	String ID2=m2.group(2);
	        	String type=m2.group(3);
	        	if(!ht_id2name.containsKey(ID1)||!ht_id2name.containsKey(ID2))
	        	{
	        		continue;
	        	}
	        	if(type.equals("maplink"))
	        	{
//	        		System.out.println("checkpoint_Methods.java_1:"+type);
	        		if((ht_id2type.get(ID1).equals("map"))&&(ht_id2type.get(ID2).equals("map")))
	        		{
	        			continue;
//	        			relationString.append(ht_id2name.get(ID1)+"\t"+"maplink"+"\t"+ht_id2name.get(ID2)+"\t"+pname+"\r\n");
	        		}
//	        			
	        		else if((ht_id2type.get(ID1).equals("map"))&&(ht_id2type.get(ID2).equals("gene")))
	        		{
//	        			System.out.println("checkpoint_Methods.java_1:"+ht_id2type.get(ID1));
	        			String reqPath="path:"+pname;
	        			String resFileName=pathdir+ht_id2name.get(ID1).split(":")[1]+".xml";
	        			String resPath=ht_id2name.get(ID1);
	        			String reqGeneName=ht_id2name.get(ID2);
	        			String cname=this.getCpdName(reader, ht_id2name);
	        			File tempFile = new File(resFileName);
	        			boolean exists = tempFile.exists();
//	        			System.out.println("checkpoint_Methods.java_1:"+reqPath+" "+resFileName+" "+resPath+" "+reqGeneName+" "+cname);
	        			if(exists)
	        				relationString.append(searchMap(reqPath,resFileName,resPath,reqGeneName,cname,"map1",ht_abst2rndnum));
	        					
	        		}
	        		else if((ht_id2type.get(ID2).equals("map"))&&(ht_id2type.get(ID1).equals("gene")))
	        		{
	        			String reqPath="path:"+pname;
	        			String resFileName=pathdir+ht_id2name.get(ID2).split(":")[1]+".xml";
	        			String resPath=ht_id2name.get(ID2);
	        			String reqGeneName=ht_id2name.get(ID1);
	        			String cname=this.getCpdName(reader, ht_id2name);
	        			File tempFile = new File(resFileName);
	        			boolean exists = tempFile.exists();
//	        			System.out.println("checkpoint_Methods.java_1:"+reqPath+" "+resFileName+" "+resPath+" "+reqGeneName+" "+cname);
	        			if(exists)
	        				relationString.append(searchMap(reqPath,resFileName,resPath,reqGeneName,cname,"map2",ht_abst2rndnum));
	        		}
	        		
	        	}

	        	
	        }
	        
	    }
		//return relationString.toString();
		return relationString.toString();
	}
	public String pathlinkfuc1(String inputfile,String pname,String pathdir) throws IOException
	{
		File rfile=new File(inputfile);
		BufferedReader reader = new BufferedReader(new FileReader(rfile));

		Hashtable<String,String> ht_id2name=new Hashtable<String,String>();    //id-rawname
		Hashtable<String,String> ht_id2type=new Hashtable<String,String>();    //id-type
		Hashtable<String,Integer> ht_abst2rndnum=new Hashtable<String,Integer>();    //inst-1, check if it is already written, for confi gene
		StringBuilder relationString=new StringBuilder();   //expandable variable, used to output results£¬ for PN, coresponding relation
		
		String tempString = null;
		while ((tempString = reader.readLine()) != null) 
	    {
	    	Pattern pattern1=Pattern.compile("<entry id=\"(.*)\".*name=\"(.*)\" type=\"([^\"]*)\"");
	        Matcher m1=pattern1.matcher(tempString);
	        Pattern pattern2=Pattern.compile("<relation entry1=\"(.*)\" entry2=\"(.*)\" type=\"(.*)\""); 
	        Matcher m2=pattern2.matcher(tempString);
	        
	        //match entry
	        if(m1.find())
	        {
	        	String id=m1.group(1);
	        	String name=m1.group(2);
	        	String type=m1.group(3);
	        	if(type.equals("gene"))
	        		metGene(ht_id2name,ht_id2type,id,name,type);
	        	else if(type.equals("compound"))
	        		metCompound(ht_id2name,ht_id2type,id,name,type);
	        	else if(type.equals("map"))
	        		metMap(ht_id2name,ht_id2type,id,name,type);
//	        	System.out.println("checkpoint_Methods.java_1:"+id);
	        }
	        
	        //match relation
	        else if(m2.find())
	        {
	        	String ID1=m2.group(1);
	        	String ID2=m2.group(2);
	        	String type=m2.group(3);
	        	if(!ht_id2name.containsKey(ID1)||!ht_id2name.containsKey(ID2))
	        	{
	        		continue;
	        	}
	        	if(type.equals("maplink"))
	        	{
//	        		System.out.println("checkpoint_Methods.java_1:"+type);
	        		if((ht_id2type.get(ID1).equals("map"))&&(ht_id2type.get(ID2).equals("map")))
	        		{
//	        			System.out.println("checkpoint_Methods.java_1:"+type);
//	        			continue;
	        			relationString.append(ht_id2name.get(ID1)+"\t"+"maplink"+"\t"+ht_id2name.get(ID2)+"\t"+pname+"\r\n");
//	        			System.out.println(ht_id2name.get(ID1)+"\t"+"maplink"+"\t"+ht_id2name.get(ID2));
	        		}
//	        			
	        		
	        	}

	        	
	        }
	        
	    }
		//return relationString.toString();
		return relationString.toString();
	}
	public void pathlinkfuc2(String inputfile,String pname,Hashtable<String,String> gene2path) throws IOException
	{
		File rfile=new File(inputfile);
		BufferedReader reader = new BufferedReader(new FileReader(rfile));

		Hashtable<String,String> ht_id2name=new Hashtable<String,String>();    //id-rawname
		Hashtable<String,String> ht_id2type=new Hashtable<String,String>();    //id-type
		Hashtable<String,Integer> ht_abst2rndnum=new Hashtable<String,Integer>();    //inst-1, check if it is already written, for confi gene
		
		String tempString = null;
		while ((tempString = reader.readLine()) != null) 
	    {
	    	Pattern pattern1=Pattern.compile("<entry id=\"(.*)\".*name=\"(.*)\" type=\"([^\"]*)\"");
	        Matcher m1=pattern1.matcher(tempString);
	        Pattern pattern2=Pattern.compile("<relation entry1=\"(.*)\" entry2=\"(.*)\" type=\"(.*)\""); 
	        Matcher m2=pattern2.matcher(tempString);
	        
	        //match entry
	        if(m1.find())
	        {
	        	String id=m1.group(1);
	        	String name=m1.group(2);
	        	String type=m1.group(3);
	        	if(type.equals("gene"))
	        	{
	        		metGene(ht_id2name,ht_id2type,id,name,type);

	        	}


//	        	System.out.println("checkpoint_Methods.java_1:"+id);
	        }
	        
	      //match relation
	        else if(m2.find())
	        {
	        	String ID1=m2.group(1);
	        	String ID2=m2.group(2);
	        	String type=m2.group(3);
	        	if(!ht_id2name.containsKey(ID1)||!ht_id2name.containsKey(ID2))
	        	{
	        		continue;
	        	}
	        	
        		String[] set1=ht_id2name.get(ID1).split("_");
        		for(int i=0;i<set1.length;i++)
        		{
        			if(gene2path.containsKey(set1[i]))
        			{
        				String ori=gene2path.get(set1[i]);
        				String[] ori_ary=ori.split("_");
        				int l=ori_ary.length;
        				if(!pname.equals(ori_ary[l-1]))
        					gene2path.put(set1[i], ori+"_"+pname);
        				
        			}
        			else
        			{
        				gene2path.put(set1[i], pname);
        				
        			}
        		}
        		
        		String[] set2=ht_id2name.get(ID2).split("_");
        		for(int i=0;i<set2.length;i++)
        		{
        			if(gene2path.containsKey(set2[i]))
        			{
        				String ori=gene2path.get(set2[i]);
        				String[] ori_ary=ori.split("_");
        				int l=ori_ary.length;
        				if(!pname.equals(ori_ary[l-1]))
        					gene2path.put(set2[i], ori+"_"+pname);
        				
        			}
        			else
        			{
        				gene2path.put(set2[i], pname);
        				
        			}
        		}


	        	
	        }
	        
	    }

	}
	public void readPathfile4MN(String inputdir,String pname) throws IOException
	{
		
	}
	public void readPathfile4PMN(String inputdir,String pname) throws IOException
	{
		
	}
}
