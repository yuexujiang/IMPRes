package Algorithm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Attribute {
	private Hashtable<String, String> kegg2id = new Hashtable<String, String>();
	public void tem_show(ArrayList<Node> end)
	{
		for(Node a : end)
		{
			String n=a.getname();
			Pattern p1=Pattern.compile(":(\\S*)", Pattern.CASE_INSENSITIVE);
			Matcher m1=p1.matcher(n);
			while(m1.find())
			{
				System.out.println(this.kegg2id.get(m1.group(1)));
			}
			
		}
	}
	public void gen_atr(Network net, FileWriter fw, ArrayList<Node> start, ArrayList<Node> end,ArrayList<Node> conf,String dir,Hashtable<String,Double> pv) throws IOException
	{
		//map should be: id->kegg
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline;
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("\t");
			if(str.length>1)
			{
			this.kegg2id.put(str[1], str[0]);
			}
		}
		
		
		for (Node v : net.getNodes().values())
		{
			String id="";
			Pattern	p1=Pattern.compile("#([^#]*)$", Pattern.CASE_INSENSITIVE);
			Matcher	m1=p1.matcher(v.getname());
				if(m1.find())
				{
					id=m1.group(1);
					if(this.kegg2id.containsKey(id))
					{
						id=this.kegg2id.get(id);
					}
				}
//			}
			
			
			
//			if((start.contains(v))&&(end.contains(v)))
//			{
//				for(Node vb: v.getContain())
//				{
//					fw.write(vb.getname()+" se "+id);
//					fw.write("\r\n");
//				}
//				
//			}
			if(start.contains(v))
			{
				for(Node vb: v.getContain())
				{
					fw.write(vb.getname()+" s "+id+" "+vb.getFc()+" "+vb.getPathway()+" "+pv.get(vb.getPathway()));
					fw.write("\r\n");
				}
					
			}
			else if(conf.contains(v))
			{
				for(Node vb: v.getContain())
				{
					fw.write(vb.getname()+" c "+id+" "+vb.getFc()+" "+vb.getPathway()+" "+pv.get(vb.getPathway()));
					fw.write("\r\n");
				}
			}
			else if(end.contains(v))
			{
				for(Node vb: v.getContain())
				{
					fw.write(vb.getname()+" e "+id+" "+vb.getFc()+" "+vb.getPathway()+" "+pv.get(vb.getPathway()));
					fw.write("\r\n");
				}
					
			}
			
			else
			{
				for(Node vb: v.getContain())
				{
					fw.write(vb.getname()+" r "+id+" "+vb.getFc()+" "+vb.getPathway()+" "+pv.get(vb.getPathway()));
					fw.write("\r\n");
				}
			}
			
		
		}
	}
	
	public void gen_atr_kegg(Network net, FileWriter fw, ArrayList<Node> start, ArrayList<Node> end, ArrayList<Node> conf) throws IOException
	{
		
		
		
		for (Node v : net.getNodes().values())
		{

			
			
			
			if((start.contains(v))&&(end.contains(v)))
			{
				for(Node vb: v.getContain())
				{
					fw.write(vb.getname()+" se");
					fw.write("\r\n");
				}
				
			}
			else if(start.contains(v))
			{
				for(Node vb: v.getContain())
				{
					fw.write(vb.getname()+" s");
					fw.write("\r\n");
				}
					
			}
			else if(end.contains(v))
			{
				for(Node vb: v.getContain())
				{
					fw.write(vb.getname()+" e");
					fw.write("\r\n");
				}
					
			}
			else if(conf.contains(v))
			{
				for(Node vb: v.getContain())
				{
					fw.write(vb.getname()+" c");
					fw.write("\r\n");
				}
			}
			else
			{
				for(Node vb: v.getContain())
				{
					fw.write(vb.getname()+" r");
					fw.write("\r\n");
				}
			}
		
		}
	}

}
