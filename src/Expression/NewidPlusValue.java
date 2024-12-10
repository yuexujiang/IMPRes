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

public class NewidPlusValue {
	public Hashtable<String, String> gmx2glyma = new Hashtable<String, String>();
	public Hashtable<String, String> glyma2gmx = new Hashtable<String, String>();
	
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
				Pattern p1=Pattern.compile("GMX:(\\d*)", Pattern.CASE_INSENSITIVE);
				Matcher m1=p1.matcher(str[0]);
				while(m1.find())
				{
					String gmx="gmx:"+m1.group(1);
					Pattern p2=Pattern.compile("GLYMA_(\\d*\\w\\d*)", Pattern.CASE_INSENSITIVE);
					Matcher m2=p2.matcher(str[1]);
					while(m2.find())
					{
						String glyma="Glyma."+m2.group(1);
						this.gmx2glyma.put(gmx, glyma);
						this.glyma2gmx.put(glyma, gmx);	
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
			int l=str[0].length();
			String glyma=str[0].substring(0, l);
			
			if(this.glyma2gmx.containsKey(glyma))
			{
				String gmx=this.glyma2gmx.get(glyma);
				sb.append(gmx);
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
	
	public static void main(String[] args) throws IOException
	{
		NewidPlusValue soy=new NewidPlusValue();
		soy.getMapping("resource/mapping.txt");	
		
		StringBuilder strb=soy.read("resource/newid_flower.txt");
		FileWriter writer=new FileWriter("resource/SOY_flower_expression.txt");
		writer.write(strb.toString());
		writer.close();
	}
	

}
