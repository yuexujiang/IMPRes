package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Get_DIPlist {

	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub
		ArrayList<String> al=new ArrayList<String>();
		FileWriter fw=new FileWriter("data/yeast/cell_wall/PPI/PPIlist.txt");
		
//		File cliVar = new File("data/Mmusc20170205CR.txt");
//		BufferedReader r1=new BufferedReader(new FileReader(cliVar));
//		String line="";
//		while((line=r1.readLine())!=null)
//		{
//			String[] str=line.split("\t");
//			String interactor1=str[0];
//			String interactor2=str[1];
//			Pattern p1=Pattern.compile("(DIP-\\d*N)", Pattern.CASE_INSENSITIVE);
//			Matcher m1=p1.matcher(interactor1);
//			if(m1.find())
//			{
//				if(!al.contains(m1.group(1)))
//				{
//					fw.write(m1.group(1)+"\r\n");
//					al.add(m1.group(1));
//				}
//				
//			}
//			Matcher m2=p1.matcher(interactor2);
//			if(m2.find())
//			{
//				if(!al.contains(m2.group(1)))
//				{
//					fw.write(m2.group(1)+"\r\n");
//					al.add(m2.group(1));
//				}
//			}
//		}
		
		
		File cliVar1 = new File("data/yeast/cell_wall/PPI/4932.protein.links.v10.5.txt");
		BufferedReader r2=new BufferedReader(new FileReader(cliVar1));
		String line2="";
		while((line2=r2.readLine())!=null)
		{
			String[] str=line2.split(" ");
			String interactor1=str[0];
			String interactor2=str[1];
			String score=str[2];
			
			if(!al.contains(interactor1))
			{
				fw.write(interactor1+"\r\n");
				al.add(interactor1);
			}
			if(!al.contains(interactor2))
			{
				fw.write(interactor2+"\r\n");
				al.add(interactor2);
			}
			
		}
		


		
		fw.close();
	}

}
