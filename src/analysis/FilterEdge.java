package analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class FilterEdge {

	public static void main(String[] args) throws IOException
	{
		// TODO Auto-generated method stub
		FileWriter writer = new FileWriter("data/lung_metastasis/filtered_net.txt");
		File pathfile=new File("data/lung_metastasis/output_net.txt");//json.json  or  TopPercent_attri.txt
		BufferedReader pathreader=new BufferedReader(new FileReader(pathfile));
		String pathline=null;
		while((pathline=pathreader.readLine())!=null)
		{
			String[] str=pathline.split("\t");
			if(!(str[1].equals("compound")))
			{
				writer.write(pathline);
				writer.write("\r\n");
			}
		}
		writer.close();
		pathreader.close();
	}

}
