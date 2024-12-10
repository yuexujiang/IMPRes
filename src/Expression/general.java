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

public class general {
	
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
		
		
//		FileWriter writer=new FileWriter("data/maize/V4/map_for_algo.txt");		
//		FileWriter writer1=new FileWriter("data/maize/V4/STRINGmap_new.txt");
		
		
		File mappingfile=new File("data/sharad_grant_mouse/uniprot_tab.txt");
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("\t");
			if((str.length==2)&&(!(str[1].equals(""))))
			{
//				Pattern	p1=Pattern.compile("// ([\\d]*) ///", Pattern.CASE_INSENSITIVE);
//				Matcher	m1=p1.matcher(str[9]);
//					if(m1.find())
					{
						String id=str[0];
						String hsa=str[1].split(";")[0];
						String kegg="abst#"+hsa;
//						String symbol=str[1].split(" /// ")[0];
						id2kegg.put(id,kegg);
//						writer.write(id+"\t"+hsa+"\r\n");
					}
							
			}
//			if((str.length==5)&&(!(str[3].equals(""))))
//			{
//					{
//						String id=str[4];
//						String hsa=str[3].split(";")[0];
//						writer.write(id+"\t"+hsa+"\r\n");
//					}
//							
//			}
//			if((!(str[2].equals("")))&&(!(str[3].equals(""))))
//			{
//				String ppi=str[3].split(";")[0];
//				String kegg=str[2].split(";")[0];
//				writer1.write(ppi+"\t"+kegg+"\r\n");
//			}
		}
		
//		writer.close();
//		writer1.close();
		
//		FileWriter writer1=new FileWriter("data/yeast/cell_wall/PPI/STRINGmap_new.txt");
//		File Stringmap=new File("data/yeast/cell_wall/PPI/STRINGfromuni.tab");     
//		BufferedReader Stringmapreader=new BufferedReader(new FileReader(Stringmap));
//		String Stringline=null;
//		while((Stringline=Stringmapreader.readLine())!=null)
//		{
//			String[] str=Stringline.split("\t");
//			if(str.length==2) // a b
//			{
//				writer.write(str[1]+"\t"+str[1]+"\r\n");
//				
//			}
//			else if(str.length==3) // a b c 
//			{
//				writer.write(str[1]+"\t"+str[1]+"\r\n");
//				writer1.write(str[1]+"\t"+str[2]+"\r\n");
//			}
//			else if(str.length==4)  
//			{
//				if(str[2].isEmpty()) // a b   d
//				{
//					writer.write(str[3].split(" ")[0]+"\t"+str[1]+"\r\n");
//				}
//				else //a b c d
//				{
//					writer.write(str[3].split(" ")[0]+"\t"+str[1]+"\r\n");
//					writer1.write(str[1]+"\t"+str[2]+"\r\n");
//				}
//			}
//		}
//		FileWriter writer2=new FileWriter("data/mouse/PPI/DIPmap_new.txt");
//		File DIPmap=new File("data/mouse/PPI/DIPfromuni.tab");
//		BufferedReader DIPmapreader=new BufferedReader(new FileReader(DIPmap));
//		String DIPline=null;
//		while((DIPline=DIPmapreader.readLine())!=null)
//		{
//			String[] str=DIPline.split("\t");
//			if(str.length==2) // a b
//			{
//				writer.write(str[1]+"\t"+str[1]+"\r\n");
//				
//			}
//			else if(str.length==3) // a b c 
//			{
//				writer.write(str[1]+"\t"+str[1]+"\r\n");
//				writer2.write(str[1]+"\t"+str[2]+"\r\n");
//			}
//			else if(str.length==4)  
//			{
//				if(str[2].isEmpty()) // a b   d
//				{
//					writer.write(str[3].split(" ")[0]+"\t"+str[1]+"\r\n");
//				}
//				else //a b c d
//				{
//					writer.write(str[3].split(" ")[0]+"\t"+str[1]+"\r\n");
//					writer2.write(str[1]+"\t"+str[2]+"\r\n");
//				}
//			}
//		}
		
//		writer1.close();
//		writer2.close();
//		writer.close();
		
		Hashtable<String, Double> kegg2value = new Hashtable<String, Double>();
		Hashtable<String, String> kegg2exp = new Hashtable<String, String>();
		
		FileWriter writer_exp=new FileWriter("data/sharad_grant_mouse/exp_sprouty2suffi_vs_del_4vs4.txt");
		
		File expfile=new File("data/sharad_grant_mouse/raw_exp_sprouty2suffi_vs_del_4vs4.txt");
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
				String[] con1=Arrays.copyOfRange(str, 1, 1+4);
				String[] con2=Arrays.copyOfRange(str, 1+4, 1+4+4);
				general x=new general();
				double[] d1=x.setExpressionfromarray(con1);
				double[] d2=x.setExpressionfromarray(con2);
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
