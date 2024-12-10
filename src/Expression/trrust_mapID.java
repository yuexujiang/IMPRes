package Expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Hashtable;

public class trrust_mapID {

	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stubHashtable<String, String> id2kegg = new Hashtable<String, String>();
		
		Hashtable<String, String> gene2pro = new Hashtable<String, String>();
//		FileWriter writer=new FileWriter("data/tarak_mouse/cell_wall/PPI/map_for_algo_add.txt");		
		File mappingfile_pro=new File("data/proteomic/gene_pro_map.txt");
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile_pro));
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
						String hsa=str[1];
						String kegg="abst#"+hsa;
//						String symbol=str[1].split(" /// ")[0];
						gene2pro.put(id,kegg);
//						writer.write(symbol+"\t"+"hsa:"+str[2].split(" /// ")[0]+"\r\n");
					}
							
			}
		}
		
		Hashtable<String, String> gene2kegg = new Hashtable<String, String>();
		File mappingfile_gene=new File("data/proteomic/gene_kegg_map.txt");
		mappingreader=new BufferedReader(new FileReader(mappingfile_gene));
		mappingline=null;
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
						String hsa=str[1];
						String kegg="abst#"+hsa;
//						String symbol=str[1].split(" /// ")[0];
						gene2kegg.put(id,kegg);
//						writer.write(symbol+"\t"+"hsa:"+str[2].split(" /// ")[0]+"\r\n");
					}
							
			}
		}
		
		
		FileWriter writer_exp=new FileWriter("data/proteomic/TF_mRNA_net.txt");
		
		File expfile=new File("data/proteomic/trrust_rawdata.human.tsv");
		BufferedReader expreader=new BufferedReader(new FileReader(expfile));
		String expline=null;
		while((expline=expreader.readLine())!=null)
		{
			String[] str=expline.split("\t");
			
			String id1=str[0];
			String id2=str[1];
			
			if((gene2pro.containsKey(id1))&&(gene2kegg.containsKey(id2)))
			{
				String pro=gene2pro.get(id1);
				String kegg=gene2kegg.get(id2);
				String rel=str[2];
				writer_exp.write(pro+"\t"+rel+"\t"+kegg+"\r\n");
			}
		

		}
		writer_exp.close();

}
}
