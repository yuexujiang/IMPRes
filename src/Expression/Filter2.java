package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Filter2 {

	Hashtable<String, Integer> dfhash=new Hashtable<>();//mmu:xxxx-->num
	Hashtable<Integer, String> num2node=new Hashtable<>();
	Hashtable<Integer, Integer> already=new Hashtable<Integer, Integer>();
	Hashtable<String, Integer> qualified=new Hashtable<>();//num-->num
	StringBuilder sb_node=new StringBuilder();
	StringBuilder sb_edge=new StringBuilder();
	
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
	
	public void readjson(String dir) throws IOException
	{
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		
		
		//{"source":1,"target":2,"relation":"expression"},
		
		int index=0;
		int highest=-1;
		while((mappingline=mappingreader.readLine())!=null)
		{
			Pattern p1=Pattern.compile("(mmu:\\d*)", Pattern.CASE_INSENSITIVE);
			Matcher m1=p1.matcher(mappingline);
			Pattern p2=Pattern.compile("\"source\":(\\d*),\"target\":(\\d*),\"relation\":(.*)}", Pattern.CASE_INSENSITIVE);
			Matcher m2=p2.matcher(mappingline);
			// \\.*relation\":(\\.*)}
			if(m1.find())
			{
				this.num2node.put(index, mappingline);
				index++;
			}
			else if(m2.find())
			{
				int num1=Integer.valueOf(m2.group(1));
				int num2=Integer.valueOf(m2.group(2));
				String relation=m2.group(3);
				System.out.println("m2group1 "+m2.group(1));
				System.out.println("m2group2 "+m2.group(2));
				System.out.println("m2group3 "+m2.group(3));

				String node1=this.num2node.get(num1);
				String node2=this.num2node.get(num2);
				Pattern pp=Pattern.compile("(mmu:\\d*)", Pattern.CASE_INSENSITIVE);
				Matcher m3=pp.matcher(node1);
				System.out.println(node1);
				Matcher m4=pp.matcher(node2);
				String bl1=null;
				String bl2=null;
				if(m3.find())
				{
					bl1=m3.group(1);
					System.out.println(bl1);
				}
				if(m4.find())
				{
					bl2=m4.group(1);
					System.out.println(bl2);
				}
				
				
				System.out.println(relation);
				
				if(this.dfhash.containsKey(bl1)||this.dfhash.containsKey(bl2))
				{
					System.out.println("a");
					if(this.already.containsKey(num1))
					{
						sb_node.append(node2+"\r\n");
						int n1=highest;
						int n2=++highest;
//						{"source":0,"target":1,"relation":"inhibition"},
						sb_edge.append("{\"source\":"+n1+",\"target\":"+n2+",\"relation\":"+relation+"},"+"\r\n");
						this.already.put(num2, 1);
					}
					else {
						sb_node.append(node1+"\r\n");
						sb_node.append(node2+"\r\n");
						int n1=++highest;
						int n2=++highest;
						sb_edge.append("{\"source\":"+n1+",\"target\":"+n2+",\"relation\":"+relation+"},"+"\r\n");
						this.already.put(num2, 1);
					}
				}
			}
		}
		mappingreader.close();
	}
	
	public static void main(String[] args) throws IOException
	{
		Filter2 ft=new Filter2();
		FileWriter fw=new FileWriter("data/mouse/dfn.json");
		
		ft.generateDFHash("data/mouse/DFall.txt");
		ft.readjson("data/mouse/json_n7.json");
		
		fw.write("{\r\n\"nodes\":[\r\n");
		fw.write(ft.sb_node.toString());
		fw.write("],\"edges\":[\r\n");
		fw.write(ft.sb_edge.toString());
		fw.write("]\r\n}\r\n");
	
	
		fw.close();
	}
	
}