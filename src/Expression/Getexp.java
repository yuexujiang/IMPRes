package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Map;

import org.apache.commons.math3.stat.inference.TTest;

public class Getexp {

	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub
		Hashtable<String, String> id2kegg = new Hashtable<String, String>();
		ArrayList<String> check=new ArrayList<String>();
		String dir="";
		String targetdir="";
		if(args.length==0)
		{
			dir="data/human_non_smoking_lung/";
			targetdir="data/human_non_smoking_lung/";
		}
		else
		{
			dir=args[0];
			targetdir=args[1];
		}
		
		File mappingfile=new File(dir+"map_for_algo.txt");
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("\t");
			if(str.length==2)
			{
				String kegg=str[1];
				String id=str[0];
				if(id2kegg.containsKey(id))
				{
					id2kegg.put(id, id2kegg.get(id)+";"+kegg);
				}
				else
				{
					id2kegg.put(id, kegg);
				}			
			}
		}
		
		Hashtable<String, Double> kegg2value = new Hashtable<String, Double>();
		Hashtable<String, String> kegg2exp = new Hashtable<String, String>();
		
		FileWriter writer_exp=new FileWriter(targetdir+"rawexp_converted.txt");
		
		File expfile=new File(targetdir+"rawexp.txt");
		BufferedReader expreader=new BufferedReader(new FileReader(expfile));
		String expline=null;
		while((expline=expreader.readLine())!=null)
		{
			String[] str=expline.split("\t");
			
			String id=str[0];
			
			if(id2kegg.containsKey(id))
			{
				String[] str2=id2kegg.get(id).split(";");
		//		String kegg=id2kegg.get(id);
				String exp=expline.split("\t", 2)[1];
				for(int i=0;i<str2.length;i++)
				{
					String kegg="abst#"+str2[i];
					kegg2exp.put(kegg, exp);
					
				}
				
				
			
			}
		}
		
		for(Map.Entry<String, String> i: kegg2exp.entrySet())
		{
			writer_exp.write(i.getKey()+"\t"+i.getValue()+"\r\n");
		}
		
		writer_exp.close();
		
		
		
		
	}

}
