package analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class KeggID2old {
	public Hashtable<String, String> gmx2glyma = new Hashtable<String, String>();
	
	public void getOldMapping(String dir) throws IOException
	{
		
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
	
				Pattern p1=Pattern.compile("(gmx:\\d*),(Glyma.*)$", Pattern.CASE_INSENSITIVE);
				Matcher m1=p1.matcher(mappingline);
				while(m1.find())
				{
						
					this.gmx2glyma.put(m1.group(1), m1.group(2));
							
				}
					
		}
		mappingreader.close();
		
	}

	public void convert(String dir) throws IOException
	{
		FileWriter writer = new FileWriter(dir+"oldID.txt");
		File pathfile=new File(dir+"KeggID.txt");//json.json  or  TopPercent_attri.txt
		BufferedReader pathreader=new BufferedReader(new FileReader(pathfile));
		String pathline=null;
		while((pathline=pathreader.readLine())!=null)
		{

			Pattern p1=Pattern.compile("(gmx:\\d*)$", Pattern.CASE_INSENSITIVE);
			Matcher m1=p1.matcher(pathline);
			if(m1.find())
			{
	        if(this.gmx2glyma.containsKey(m1.group(1)))
	        {
	        	
	        	writer.write(this.gmx2glyma.get(m1.group(1)));
	        	writer.write("\r\n");

	        }
			}
		}
		writer.close();
	}
	
	public static void main(String[] args) throws IOException
	{
		KeggID2old a=new KeggID2old();
		a.getOldMapping("resource/mapping_ulti.txt");
		a.convert("resource/result_gmx_RHheat_nonode/");
	}
}
