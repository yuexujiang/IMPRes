package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Filter {
	Hashtable<String, Integer> dfhash=new Hashtable<>();//mmu:xxxx-->num
	Hashtable<String, Integer> qualified=new Hashtable<>();//num-->num
	
	public void generateDFHash(String dir) throws IOException
	{
		
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
				Pattern p1=Pattern.compile("(mmu:\\d*)", Pattern.CASE_INSENSITIVE);
				Matcher m1=p1.matcher(mappingline);
				while(m1.find())
				{
					String id=m1.group(1);
					this.dfhash.put(id, 1);
				}
					
		}
		mappingreader.close();
	}
	
	public String readjson(String dir) throws IOException
	{
		StringBuilder sb=new StringBuilder();
		File f=new File(dir);
		BufferedReader bfreader=new BufferedReader(new FileReader(f));
		String line=null;
		int index=0;
//		{
//			"nodes":[
//		sb.append("{\r\n\"nodes\":[\r\n");
		while((line=bfreader.readLine())!=null)
		{
			Pattern p1=Pattern.compile("(mmu:\\d*)", Pattern.CASE_INSENSITIVE);
			Matcher m1=p1.matcher(line);
			Pattern p2=Pattern.compile("\"source\":(\\d*),\"target\":(\\d*)", Pattern.CASE_INSENSITIVE);
			Matcher m2=p2.matcher(line);
			if(m1.find())
			{
				String id=m1.group(1);
				if(this.dfhash.containsKey(id))
				{
					String ind=String.valueOf(index);
					this.qualified.put(ind, 1);
					
					sb.append(line+"\r\n");
				}
				index++;
			}
			else if(m2.find())
			{
				String key1=m2.group(1);
				String key2=m2.group(2);
				if(this.qualified.containsKey(key1)&&this.qualified.containsKey(key2))
				{
					sb.append(line+"\r\n");
				}
			}
			else {
				sb.append(line+"\r\n");
			}
		}
		return sb.toString();
	}
			
		

	public static void main(String[] args) throws IOException
	{
		FileWriter fw=new FileWriter("data/mouse/dfwtjson.json");
		Filter ft=new Filter();
		ft.generateDFHash("data/mouse/DFall.txt");
		String content=ft.readjson("data/mouse/json_wt4.json");
		fw.write(content);
		fw.close();
	}
}
