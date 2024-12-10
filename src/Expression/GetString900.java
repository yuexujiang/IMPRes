package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GetString900 {

	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub
		FileWriter fw=new FileWriter("data/Ron_soybean_flower/PPI/String950.txt");
		
		File StringNET=new File("data/Ron_soybean_flower/PPI/3847.protein.links.v11.5.txt");
		BufferedReader reader = new BufferedReader(new FileReader(StringNET));
		reader.readLine();
		String line="";
		int check=0;
		while ((line = reader.readLine()) != null)
		{
			check++;
			String[] str=line.split(" ");
			if(Integer.valueOf(str[2])<950)
				continue;
			else
			{
				fw.write(line+"\r\n");
			}
			if(check%10000==0)
			{
				System.out.println(check);
			}
		}
		
		fw.close();

	}

}
