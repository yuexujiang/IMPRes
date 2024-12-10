package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.inference.TTest;

public class general_convert {
	
	public double[] setExpressionfromarray(String[] db)
	{
		int s=db.length;
		double[] a=new double[s];
			for(int j=0;j<s;j++)
			{
				a[j]=Double.valueOf(db[j]);
			}
		return a;
	
	}

	
	public static void main(String[] args) throws IOException{
	
		
		Hashtable<String, String> id2kegg = new Hashtable<String, String>();
		
		
		File mappingfile=new File("data/Hekmat_human/uniprot_mapping.tab");
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			System.out.println(mappingline);
			
			String[] str=mappingline.split("\t");
			System.out.println(str.length);
//			if((str.length==4)&&((str[2].substring(0, 3).equals("hsa"))))
			if((str.length==4)&&((str[2]!="")))
			{
					{
						String id=str[1];
						String hsa=str[2].split(";")[0];
						String kegg="abst#"+hsa;
//						String symbol=str[1].split(" /// ")[0];
						id2kegg.put(id,kegg);
//						writer.write(id+"\t"+hsa+"\r\n");
					}
							
			}

		}
		
		
		Hashtable<String, Double> kegg2value = new Hashtable<String, Double>();
		Hashtable<String, String> kegg2exp = new Hashtable<String, String>();
		
		FileWriter writer_exp=new FileWriter("data/Hekmat_human/exp.txt");
		
		File expfile=new File("data/Hekmat_human/raw_exp.txt");
		BufferedReader expreader=new BufferedReader(new FileReader(expfile));
		String expline=null;
		double s=1.0;
		while((expline=expreader.readLine())!=null)
		{
			String[] str=expline.split("\t");
			
			String id=str[0];
			
			if(id2kegg.containsKey(id))
			{
				String kegg=id2kegg.get(id);
				String exp=expline.split("\t", 2)[1];
				String[] con1=Arrays.copyOfRange(str, 1, 1+3);
				String[] con2=Arrays.copyOfRange(str, 1+3, 1+3+3);
				general x=new general();
				double[] d1=x.setExpressionfromarray(con1);
				double[] d2=x.setExpressionfromarray(con2);
				double sumd1=0.0;
				for(double i:d1)
				{
					sumd1=sumd1+i;
				}
				double sumd2=0.0;
				for(double i:d2)
				{
					sumd2=sumd2+i;
				}
				if((sumd1==0.0)||(sumd2==0.0))
					continue;
				
				s=new TTest().tTest(d1, d2);
				if(kegg2value.containsKey(kegg))
				{
					if(s<kegg2value.get(kegg))
					{
						kegg2value.put(kegg, s);
						kegg2exp.put(kegg, exp);
					}
					
				}
				else
				{
					kegg2value.put(kegg, s);
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
