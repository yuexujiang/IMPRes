package kegg;
import java.io.*;
import java.util.*;
import java.util.regex.*;

public class Convert 
{
	public static void main(String[] args) throws IOException
	{
		String indir="data/poplar/net/";	//input directory
//		String code="gmx";	//species code	
//		ArrayList<Node> finalList=new ArrayList<Node>();//store the final nodes information
		
		//specify the output file addresses and initial them
		//four output files for network, start, end and confidence respectively
		String output_net_dir=indir+"output_net.txt";
		String output_start_dir=indir+"output_start.txt";
		String output_end_dir=indir+"output_end.txt";
		String output_maplink_dir=indir+"output_map.txt";
		String output_pathlink_dir=indir+"output_path.txt";
//		String output_conf_dir=indir+"output_conf.txt";
		FileWriter PMNwriter = new FileWriter(output_net_dir);
		FileWriter PMNstart=new FileWriter(output_start_dir);
		FileWriter PMNend=new FileWriter(output_end_dir);
		FileWriter maplink=new FileWriter(output_maplink_dir);
		FileWriter pathlink=new FileWriter(output_pathlink_dir);
		
		FileWriter MNwriter = new FileWriter(indir+"output_MNnet.txt");
		FileWriter MNstart=new FileWriter(indir+"output_MNstart.txt");
		FileWriter MNend=new FileWriter(indir+"output_MNend.txt");
//		FileWriter conf=new FileWriter(output_conf_dir);
		Hashtable<String,String> ht_gene2path=new Hashtable<String,String>();
		
		Methods metho=new Methods();
		//find all the pathways, print how many pathways we have
		File[] tempList = metho.getPathfiles(indir+"kgml");
		int N=tempList.length;
		System.out.println("pathway number N is��"+N);
		
		
		for(int z=0;z<tempList.length;z++)
		{
			//get pathway name, and print the ongoing pathway...
			String pathname="";
			Pattern p=Pattern.compile("(\\w*)\\.xml");
	        Matcher m=p.matcher(tempList[z].toString());
	        if(m.find())
	        {
	        	pathname=m.group(1).toString();
	        }			
			System.out.println("parsing..."+pathname);
			
			//call readPathfile function
			String input_file=tempList[z].toString();
			metho.readPathfile4PPN(input_file, pathname);//read pathway file one by one
			
			String path=metho.relation;			
			PMNwriter.write(path);//record current relation information
			
			path=metho.reaction;
			MNwriter.write(path);
	//		finalList.addAll(metho.al);//record current pathway nodes' information
			
	//		metho.al.clear();//clear arraylist
			metho.relation="";//clear String
			metho.reaction="";
			
			String pathdir=indir+"kgml/";
			String maplink_str=metho.readmaplink(input_file,pathname,pathdir);
			maplink.write(maplink_str);
			
			String pathlink_str=metho.pathlinkfuc1(input_file,pathname,pathdir);
			pathlink.write(pathlink_str);
			metho.pathlinkfuc2(input_file,pathname,ht_gene2path);
			
		}
		for(String g : ht_gene2path.keySet())
		{
			System.out.println(g);
			System.out.println(ht_gene2path.get(g));
			String[] set1=ht_gene2path.get(g).split("_");
			int n=set1.length;
			for(int i=0;i<n-1;i++)
			{
				for(int j=i+1;j<n;j++)
				{
					pathlink.write(set1[i]+"\tmaplink\t"+set1[j]+"\r\n");
					pathlink.write(set1[j]+"\tmaplink\t"+set1[i]+"\r\n");
				}
				
			}
		}
		
		
		metho.writeHash(metho.htIndegree, PMNstart);
		metho.writeHash(metho.htOutdegree, PMNend);
		metho.writeHash(metho.htIndegree_MN, MNstart);
		metho.writeHash(metho.htOutdegree_MN, MNend);
		PMNwriter.close();
		PMNstart.close();
		PMNend.close();
		maplink.close();
		pathlink.close();
		
		MNwriter.close();
		MNstart.close();
		MNend.close();
	}
}
