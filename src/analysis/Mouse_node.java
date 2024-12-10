package analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class Mouse_node {

	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub
		Hashtable<String, String> id2func=new Hashtable();
		Hashtable<String, String> path2pname=new Hashtable();
		
		File pathmapping=new File("data/mouse/Allgenes.txt");
		BufferedReader mappingreader=new BufferedReader(new FileReader(pathmapping));
		String mappingline;
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("\t");
			if(str.length==2)
			{
				
				String id=str[0];
				String func=str[1];
				id2func.put(id, func);
			}
		}
		
		pathmapping=new File("data/mouse/mouse_pathway map.txt");
		mappingreader=new BufferedReader(new FileReader(pathmapping));
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("  ");
			if(str.length==2)
			{
				
				String path="mmu"+str[0];
				String pname=str[1];
				path2pname.put(path, pname);
			}
		}
		
		FileWriter fw=new FileWriter("data/mouse/New folder/lps_info.txt");
		
		File inputfile=new File("data/mouse/New folder/lps_node.txt");
		BufferedReader inputreader=new BufferedReader(new FileReader(inputfile));
		String inputline;
		while((inputline=inputreader.readLine())!=null)
		{
			String[] str=inputline.split("\t");
			if(str.length==8)
			{
				String id=str[2].toUpperCase();
				String fd=str[3];
				String path=str[4];
				fw.write(id+"\t"+fd+"\t"+path+"\t"+path2pname.get(path)+"\t"+id2func.get(id)+"\r\n");
			}
			
		}
		
		fw.close();
		
	}

}
