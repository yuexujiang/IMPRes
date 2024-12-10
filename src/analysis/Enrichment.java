package analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Enrichment 
{
	public Hashtable<String, Integer> ht_path=new Hashtable();
	public int Nodenum=0;
	public int hitnum=0;
	
	public void nodeenrich(String dir) throws IOException
	{
		FileWriter writer = new FileWriter(dir+"keggID.txt");
		File pathfile=new File(dir+"json.json");//json.json  or  TopPercent_attri.txt
		BufferedReader pathreader=new BufferedReader(new FileReader(pathfile));
		String pathline=null;
		while((pathline=pathreader.readLine())!=null)
		{
			
			Pattern pathpattern=Pattern.compile("(gmx:\\d+).*(gmx\\d+)\"");
	        Matcher pathm=pathpattern.matcher(pathline);
	        if(pathm.find())
	        {
	        	String keggid=pathm.group(1);
	        	String pathid=pathm.group(2);
	        	writer.write(keggid);
	        	writer.write("\r\n");
	        	if(this.ht_path.containsKey(pathid))
	        	{
	        		this.hitnum++;
	        	}
	        	Nodenum++;
	        }
		}
		writer.close();
	}
	
	public void pathEnrich(String dir)
	{
		
	}
	
	public Hashtable<String, Integer> gethtpath(String dir,String code) throws IOException
	{
		File pathfile=new File(dir);
		BufferedReader pathreader=new BufferedReader(new FileReader(pathfile));
		String pathline=null;
		while((pathline=pathreader.readLine())!=null)
		{
			
			Pattern pathpattern=Pattern.compile("^(\\d+)");
	        Matcher pathm=pathpattern.matcher(pathline);
	        if(pathm.find())
	        {
	        	this.ht_path.put(code+pathm.group(1), 1);
	        }
		}
		return this.ht_path;
	}
	
	public static void main(String[] args) throws IOException
	{
		Enrichment enri=new Enrichment();
		enri.gethtpath("resource/result_gmx_RHheat_nonode/path.txt", "gmx");
		System.out.println(enri.ht_path.size());
		enri.nodeenrich("resource/result_gmx_RHheat_nonode/");
		System.out.println(enri.Nodenum+" "+enri.hitnum);
	}
}
