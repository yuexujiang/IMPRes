package analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Drawpath {
	private Hashtable<String, String> ht_allpath=new Hashtable();
	
	public Hashtable<String, String> gethtallpath(String dir,String code)
	{
		try
		{
			
			
			File pathfile=new File(dir);
			BufferedReader pathreader=new BufferedReader(new FileReader(pathfile));
			String pathline=null;
			while((pathline=pathreader.readLine())!=null)
			{
				
				Pattern pathpattern=Pattern.compile("(\\d{5})\\s+(.*)$");
    	        Matcher pathm=pathpattern.matcher(pathline);
    	        if(pathm.find())
    	        {
    	     //   	System.out.println(pathline);
    	        	this.ht_allpath.put(code+pathm.group(1), pathm.group(2));
    	        }
			}
		//	System.out.println(this.ht_path.size());
			return this.ht_allpath;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return this.ht_allpath;
		}
		
	}
	
	public void getscore(String dir)
	{
		int through=0;
		int relevant=0;
		try
		{
			File pathfile=new File(dir);
			BufferedReader pathreader=new BufferedReader(new FileReader(pathfile));
			String pathline=null;
			FileWriter writer = new FileWriter(dir+"output3.txt");
			
			pathline=pathreader.readLine();
			Pattern pathpattern=Pattern.compile("(hsa\\d+).*(hsa\\d+)");
	        Matcher pathm=pathpattern.matcher(pathline);
	        
	        String last="aa";
	        int lastn=1;
	        if(pathm.find())
	        {
	        	System.out.println("1");
	        	String path1=pathm.group(1);
	        	String path2=pathm.group(2);
	//        	writer.write(path1+"\r\n"+path2+"\r\n");
	        	
	        	   last=path1;
	        	   if(path2.equals(last))
	        		   lastn=lastn+1;
	        	   else
	        	   {
	        		   writer.write(last+"\t"+last+": "+ht_allpath.get(last)+"\t"+lastn+"\r\n");
	        		   lastn=1;
	        		   last=path2;
	        	   }
	
	        }
			
	
			while((pathline=pathreader.readLine())!=null)
			{
				System.out.println("2");
				pathpattern=Pattern.compile("(hsa\\d+).*(hsa\\d+)");
    	        pathm=pathpattern.matcher(pathline);
    	        if(pathm.find())
    	        {
    	        	String path1=pathm.group(1);
    	        	String path2=pathm.group(2);
   // 	        	if(!path1.equals(path2))
   // 	        	writer.write(path2+"\r\n");
    	        	
    	        	last=path1;
 	        	   if(path2.equals(last))
 	        		   lastn=lastn+1;
 	        	   else
 	        	   {
 	        		   writer.write(last+"\t"+last+": "+ht_allpath.get(last)+"\t"+lastn+"\r\n");
 	        		   lastn=1;
 	        		   last=path2;
 	        	   }

    	        }
			}
			writer.write(last+"\t"+last+": "+ht_allpath.get(last)+"\t"+lastn+"\r\n");
	        System.out.println(through+" "+relevant);
	        writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
			
		}
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Drawpath a=new Drawpath();
		Hashtable<String, String> ha=new Hashtable();
    	ha=a.gethtallpath("resource/allpathes_hsa.txt","hsa");
    	
    	System.out.println(ha.get("hsa00030"));
    	a.getscore("resource/best_path_hsa.txt");

	}

}
