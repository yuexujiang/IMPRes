package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class B_ALL {
	public Hashtable<String, String> kegg2id = new Hashtable<String, String>();
	public Hashtable<String, String> id2kegg = new Hashtable<String, String>();
	
	public void getMapping(String dir) throws IOException
	{
		
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("\t");
			if(str.length==2)
			{
				Pattern p1=Pattern.compile("(^.*$)", Pattern.CASE_INSENSITIVE);
				Matcher m1=p1.matcher(str[0]);
				while(m1.find())
				{
					String id=m1.group(1);
					Pattern p2=Pattern.compile("(^\\d*)_at", Pattern.CASE_INSENSITIVE);
					Matcher m2=p2.matcher(str[1]);
					while(m2.find())
					{
						String kegg="abst#hsa:"+m2.group(1);
						this.kegg2id.put(kegg, id);
						this.id2kegg.put(id, kegg);	
					}
					
				}
			}
			
		}
		mappingreader.close();
		
	}

	public StringBuilder read(String dir) throws IOException
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
			
	
	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub
		B_ALL h=new B_ALL();
		h.getMapping("data/B_ALL/mapping.txt");	
		StringBuilder ss=h.read("data/B_ALL/name_exp.txt");
		FileWriter writer=new FileWriter("data/B_ALL/kegg_exp.txt");
		writer.write(ss.toString());
		writer.close();
	}

}
