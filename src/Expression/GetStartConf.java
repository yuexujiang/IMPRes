package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

public class GetStartConf {
	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub
		Hashtable<String, String> id2kegg = new Hashtable<String, String>();
		ArrayList<String> check=new ArrayList<String>();
		GetStartConf mn=new GetStartConf();
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
		
		
		
		FileWriter writer=new FileWriter(targetdir+"rawlist_converted.txt");
		
		
		File genefile=new File(targetdir+"rawlist.txt");
		BufferedReader genereader=new BufferedReader(new FileReader(genefile));
		String geneline=null;
		while((mappingline=genereader.readLine())!=null)
		{
			String id=mappingline;
		//	System.out.println(id);
			if(id2kegg.containsKey(id))
			{
		//		System.out.println(id);
		//		System.out.println(id2kegg.get(id));
				String[] str=id2kegg.get(id).split(";");
				
				for(int i=0;i<str.length;i++)
				{
					writer.write("abst#"+str[i]+"\r\n");
				}
			}
		}
		
		writer.close();
		
	}

}
