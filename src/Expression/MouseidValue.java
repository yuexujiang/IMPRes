package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MouseidValue {
	public Hashtable<String, String> kegg2id = new Hashtable<String, String>();
	public Hashtable<String, String> id2kegg = new Hashtable<String, String>();
	public Hashtable<String, String> name2id = new Hashtable<String, String>();
	public Hashtable<String, String> id2name = new Hashtable<String, String>();
	
	public void getGnsGeneMapping(String dir) throws IOException
	{
		
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("\t");
			if(str.length==2)
			{
				Pattern p1=Pattern.compile("(ENSMUSG\\d*)", Pattern.CASE_INSENSITIVE);
				Matcher m1=p1.matcher(str[0]);
				while(m1.find())
				{
					String id=m1.group(1);
					Pattern p2=Pattern.compile("(\\S*)", Pattern.CASE_INSENSITIVE);
					Matcher m2=p2.matcher(str[1]);
					while(m2.find())
					{
						String name=m2.group(1);
						this.name2id.put(name, id);
						this.id2name.put(id, name);	
					}
					
				}
			}
			
		}
		mappingreader.close();
		
	}
	
	public void getMapping(String dir) throws IOException
	{
		
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("\t");
			if(str.length==3)
			{
				Pattern p1=Pattern.compile("(ENSMUSG\\d*)", Pattern.CASE_INSENSITIVE);
				Matcher m1=p1.matcher(str[1]);
				while(m1.find())
				{
					String id=m1.group(1);
					Pattern p2=Pattern.compile("(mmu:\\d*)", Pattern.CASE_INSENSITIVE);
					Matcher m2=p2.matcher(str[2]);
					while(m2.find())
					{
						String kegg="abst#"+m2.group(1);
						this.kegg2id.put(kegg, id);
						this.id2kegg.put(id, kegg);	
					}
					
				}
			}
			
		}
		mappingreader.close();
		
	}
	
	public StringBuilder readid(String dir) throws IOException
	{
		StringBuilder sb=new StringBuilder();
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			
			String[] str=mappingline.split("\t");
			
			String id=str[0];
			
			if(this.id2kegg.containsKey(id))
			{
				String kegg=this.id2kegg.get(id);
				sb.append(kegg);
				for(int i=1;i<str.length;i++)
				{
					sb.append("\t"+str[i]);
				}
				sb.append("\r\n");
			}
		}
		mappingreader.close();
		return sb;
	}
	
	public StringBuilder readname(String dir) throws IOException
	{
		StringBuilder sb=new StringBuilder();
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			
			String[] str=mappingline.split("\t");
			
			String name=str[0];
			
			if(this.name2id.containsKey(name))
			{
				String id=this.name2id.get(name);
				if(this.id2kegg.containsKey(id))
				{
					String kegg=this.id2kegg.get(id);
					sb.append(kegg);
					for(int i=1;i<str.length;i++)
					{
						sb.append("\t"+str[i]);
					}
					sb.append("\r\n");
				}		
			}
		}
		mappingreader.close();
		return sb;
	}
	
	public static void main(String[] args) throws IOException
	{
		MouseidValue soy=new MouseidValue();
		soy.getMapping("data/mouse/mousemapping.tab");	
		soy.getGnsGeneMapping("data/mouse/mouseEnsGNameMapping.txt");
		
		StringBuilder strb=soy.readid("data/mouse/mouseid_value_n.txt");		
		FileWriter writer=new FileWriter("data/mouse/mouse_expression_n2.txt");
		writer.write(strb.toString());
		writer.close();
		
		StringBuilder strb2=soy.readname("data/mouse/DFgenes.txt");
		FileWriter writer2=new FileWriter("data/mouse/DFstart.txt");
		writer2.write(strb2.toString());
		writer2.close();
	}
	

}
