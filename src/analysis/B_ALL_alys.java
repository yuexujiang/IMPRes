package analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import Expression.B_ALL;

public class B_ALL_alys {

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
						String kegg=m2.group(1);
						this.kegg2id.put(kegg, id);
						this.id2kegg.put(id, kegg);	
					}
					
				}
			}
			
		}
		mappingreader.close();
		
	}
	
	public void read(String dir) throws IOException
	{
		StringBuilder sb=new StringBuilder();
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
	//		System.out.println("1");
			Pattern p2=Pattern.compile("hsa:(\\d*)", Pattern.CASE_INSENSITIVE);
			Matcher m2=p2.matcher(mappingline);
			while(m2.find())
			{
			//	System.out.println("1");
				System.out.println(this.kegg2id.get(m2.group(1)));
			}
			
			
			
		}
		mappingreader.close();
	
	}
	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub
		B_ALL_alys h=new B_ALL_alys();
		h.getMapping("data/B_ALL/mapping.txt");	
		h.read("data/B_ALL/output3.txt");
//		FileWriter writer=new FileWriter("data/B_ALL/kegg_exp.txt");
//		writer.write(ss.toString());
//		writer.close();
	}
}
