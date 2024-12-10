package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class HumLungAccession {

	public StringBuilder read(String dir) throws IOException
	{
		StringBuilder sb=new StringBuilder();
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			
			String[] str=mappingline.split("\t");
			
			Pattern p1=Pattern.compile("(^\\d*)_at", Pattern.CASE_INSENSITIVE);
			Matcher m1=p1.matcher(str[2]);
			while(m1.find())
			{
				sb.append("hsa:"+m1.group(1));
				for(int i=3;i<str.length;i++)
				{
					sb.append("\t"+str[i]);
				}
				sb.append("\r\n");
			}
			
		}
		return sb;
	}
			
	
	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub
		HumLungAccession h=new HumLungAccession();
		StringBuilder ss=h.read("resource/HumLungAccession.txt");
		FileWriter writer=new FileWriter("resource/HumLungAccession_expression.txt");
		writer.write(ss.toString());
		writer.close();
	}

}
