package Algorithm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;


//run sample: data/nsclc/hcc827vsh820/ unisam_self_tfun1 1 1 output_net.txt conf.txt start.txt exp.txt result_TlistTF2.txt end.txt

public class Main
{
	private String inputfolder="data/";
	private String dataType="two_conditions"; //2c_self_tfun1;2c_list_tfun1;time_list_tfun2;no_list_tfunstep;unisam_self_tfununi;dfgiven_tfun1;
	private int Num_sample1=10;
	private int Num_sample2=10;
	private String netname="output_net.txt";
	private String index_ifconf="no";
	private String start_gene="output_start.txt";
	private String index_ifexp="no";
	private String outname="result.json";
	private int Num_target=50;
	private String index_ifend="no";
	private String index_ifppi="no";
	private Hashtable<String, Double> pv=new Hashtable<String, Double>();
	private int alpha=1;
	private String TF_mRNAfile;
	private String pro_exp_location;
	private String metab_exp;
	private int num_sample1_metab=3;
	private int num_sample2_metab=3;
	
	public static void main(String args[]) throws IOException, REngineException, REXPMismatchException
	{
		Main mn = new Main();
		IO fio = new IO();
		
//		Hashtable<String, String> ctr=new Hashtable<String, String>();
//		if(args.length==0)
//		{
//			return;
//		}
//		else if(args[0].equals("-help"))
//		{
//			System.out.println("Show help");
//			return;
//		}
//		else
//		{
//			for(int ai=1;ai<args.length;ai++)
//			{
//				if(ai%2==1)
//				{
//					ctr.put(args[ai], args[ai+1]);
//					ai++;
//				}
//			}
//			if(!ctr.containsKey("-dir"))
//			{
//				System.out.println("Please specify workspace");
//				return;
//			}
//			if(!ctr.containsKey("-type"))
//			{
//				System.out.println("Please specify datatype");
//				return;
//			}
//			if(!ctr.containsKey("-net"))
//			{
//				System.out.println("Please specify kegg network");
//				return;
//			}
//			if(!ctr.containsKey("-seed"))
//			{
//				System.out.println("Please specify seed genes");
//				return;
//			}
//			if(!ctr.containsKey("-out"))
//			{
//				System.out.println("Please specify output file");
//				return;
//			}
//			
//			if(ctr.get("-type").equals("two_conditions"))
//			{
//				if(!ctr.containsKey("-n1"))
//				{
//					System.out.println("Please specify sample numbers of condition 1");
//					return;
//				}
//				if(!ctr.containsKey("-n2"))
//				{
//					System.out.println("Please specify sample numbers of condition 2");
//					return;
//				}
//				
//			}
//		}
		
		
	
		//input directory
		mn.inputfolder=args[0];
		//data type
		mn.dataType=args[1];
		//num of samples in each condition
		mn.Num_sample1=Integer.valueOf(args[2]);
		mn.Num_sample2=Integer.valueOf(args[3]);
		//network file name
		mn.netname=args[4];
		mn.index_ifend=args[10];
		mn.index_ifppi=args[11];
		mn.start_gene=args[6];
		mn.index_ifconf=args[5];
		mn.index_ifexp=args[7];
		mn.outname=args[8];
		mn.Num_target=Integer.valueOf(args[9]);
	//	mn.TF_mRNAfile=args[12];
		mn.pro_exp_location=args[12];
		mn.metab_exp=args[13];
		mn.num_sample1_metab=Integer.valueOf(args[14]);
		mn.num_sample2_metab=Integer.valueOf(args[15]);
		
		Network inNetwork=new Network();
		if(mn.dataType.equals("2c_self_tfun1"))
		{
			inNetwork=fio.readNetworkfromFile(mn.inputfolder+mn.netname, mn.Num_sample1,mn.Num_sample2);
		}
		if(mn.dataType.equals("unisam_self_tfununi"))
		{
			inNetwork=fio.readNetworkfromFile(mn.inputfolder+mn.netname, mn.Num_sample1,mn.Num_sample2);
		}
		if(mn.dataType.equals("dfgiven_tfun1"))
		{
			inNetwork=fio.readNetworkfromFile(mn.inputfolder+mn.netname, mn.Num_sample1,mn.Num_sample2);
		}
		else if(mn.dataType.equals("time_list_tfun2"))
		{
			inNetwork=fio.readNetworkfromFile(mn.inputfolder+mn.netname, mn.Num_sample1*mn.Num_sample2,mn.Num_sample1*mn.Num_sample2);
		}
		else if(mn.dataType.equals("noexp_tfun3"))
		{
			inNetwork=fio.readNetworkfromFile(mn.inputfolder+mn.netname, mn.Num_sample1,mn.Num_sample2);
		}
		else if(mn.dataType.equals("inte_pro"))
		{
			inNetwork=fio.readNetworkfromFile(mn.inputfolder+mn.netname, mn.Num_sample1,mn.Num_sample2);
			
		}
		else if(mn.dataType.equals("inte_metab"))
		{
			inNetwork=fio.readNetworkfromFile(mn.inputfolder+mn.netname, mn.Num_sample1,mn.Num_sample2);
//			unmark ID.java, redNetworkfromFile: edgetype.equals("compound")
//			unmark Presentation.java, Present: else if(interab.getType().equals("compound"))
		}
		
		if(!mn.index_ifppi.equals("no"))
		{
			ArrayList<String> al=fio.getexistingEdgelist(inNetwork);
			if(mn.dataType.equals("2c_self_tfun1"))
			{
				inNetwork=fio.readPPIfromFile(al, inNetwork, mn.inputfolder+"PPI/String950.txt", mn.inputfolder+"PPI/STRINGmap_new.txt", mn.Num_sample1,mn.Num_sample2);
			}
			if(mn.dataType.equals("unisam_self_tfununi"))
			{
				inNetwork=fio.readPPIfromFile(al, inNetwork, mn.inputfolder+"PPI/String950.txt", mn.inputfolder+"PPI/STRINGmap_new.txt", mn.Num_sample1,mn.Num_sample2);
			}
			if(mn.dataType.equals("dfgiven_tfun1"))
			{
				inNetwork=fio.readPPIfromFile(al, inNetwork, mn.inputfolder+"PPI/String950.txt", mn.inputfolder+"PPI/STRINGmap_new.txt", mn.Num_sample1,mn.Num_sample2);
			}
			else if(mn.dataType.equals("time_list_tfun2"))
			{
				inNetwork=fio.readPPIfromFile(al, inNetwork, mn.inputfolder+"PPI/String950.txt", mn.inputfolder+"PPI/STRINGmap_new.txt", mn.Num_sample1*mn.Num_sample2,mn.Num_sample1*mn.Num_sample2);
			}
			else if(mn.dataType.equals("noexp_tfun3"))
			{
				inNetwork=fio.readPPIfromFile(al, inNetwork, mn.inputfolder+"PPI/String950.txt", mn.inputfolder+"PPI/STRINGmap_new.txt", mn.Num_sample1,mn.Num_sample2);
			}
			else if(mn.dataType.equals("inte_pro"))
			{
				inNetwork=fio.readPPI4Inte(al, inNetwork, mn.inputfolder+"PPI/String950.txt", mn.Num_sample1,mn.Num_sample2);
				inNetwork=fio.readTFmRNA4Inte(al, inNetwork, mn.inputfolder+"TF_mRNA_net.txt", mn.Num_sample1,mn.Num_sample2);
			}
			
//			inNetwork=fio.readPPIfromFile(al, inNetwork, mn.inputfolder+"PPI/DIP_net.txt", mn.inputfolder+"PPI/DIP_map.tab", mn.inputfolder+"PPI/String900plusBADC500.txt", mn.inputfolder+"PPI/String_map.tab", mn.Num_sample1,mn.Num_sample2);
//			inNetwork=fio.readPPIfromFile(al, inNetwork, mn.inputfolder+"PPI/DIP_net.txt", mn.inputfolder+"PPI/DIPmap_new.txt", mn.inputfolder+"PPI/String900.txt", mn.inputfolder+"PPI/STRINGmap_new.txt", mn.Num_sample1,mn.Num_sample2);
		}
		
		//permutation
//		int pj=1;
//		{
//			
//			FileWriter per = new FileWriter(mn.inputfolder+"PERscore2.txt");
//			
//			ArrayList<Node> target=new ArrayList<Node>();
//	//		target.add(inNetwork.getByName("abst#sce:YLR332W"));
//	//		target.add(inNetwork.getByName("abst#sce:YOR008C"));
//			target.add(inNetwork.getByName("abst#sce:YPR165W"));
//			target.add(inNetwork.getByName("abst#sce:YBL105C"));
//			target.add(inNetwork.getByName("abst#sce:YLR371W"));
//	//		target.add(inNetwork.getByName("abst#sce:YBL0807"));
//			target.add(inNetwork.getByName("abst#sce:YJL095W"));
//	//		target.add(inNetwork.getByName("abst#sce:YOR231W"));
//			target.add(inNetwork.getByName("abst#sce:YPL140C"));
//			target.add(inNetwork.getByName("abst#sce:YHR030C"));
//			target.add(inNetwork.getByName("abst#sce:YPL089C"));
//	//		target.add(inNetwork.getByName("abst#sce:YER111C"));
//	//		target.add(inNetwork.getByName("abst#sce:YLR182W"));
//			
//			for(int pi=1;pi<=100000;pi++)
//			{
//				ArrayList<Node> blacklist=new ArrayList<Node>();
//				Random ran = new Random();
//				int dis=0;
//				int score=0;
//				Node Node_current=inNetwork.getByName("abst#sce:YOR008C");
//				blacklist.add(Node_current);
//				while(dis<7)
//				{
//	//				System.out.println(dis);
//					int l=Node_current.getDown().size();
//					int x = ran.nextInt(l);
//					
//							int che=0;
//							for(Node k:Node_current.getDown())
//							{
//								if(blacklist.contains(k))
//								{
//									che++;
//								}
//							}
//							if(che==l)
//								break;
//					if(blacklist.contains(Node_current.getDown().get(x)))
//					{
//						continue;
//					}
//					Node nd_next=Node_current.getDown().get(x);
//					blacklist.add(nd_next);
//					String endnodes = Node_current.getname() + nd_next.getname();
//					Interaction inter = inNetwork.getByendnames(endnodes);
//					if(inter.getType().equals("type"))
//					{
//						Node_current=nd_next;
//					}
//					else
//					{
//						Node_current=nd_next;
//						dis++;
//						if(target.contains(nd_next.getBelong()))
//						{
//							score++;
//						}
//					}
//					
//				}
//				per.write(score+"\r\n");
//				System.out.println(pi);
//			}
//			per.close();
//			if(pj==1)
//			{
//				return;
//			}
//			
//		}
		
		
		
		
		
		
		///////////////////////////////////dream part////////////////////////////////////////
		//add nodes into the network
//		File mappingfile=new File(mn.inputfolder+"map_bt20.txt");
//		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
//		String mappingline;
//		while((mappingline=mappingreader.readLine())!=null)
//		{
//			String[] str=mappingline.split("\t");
//			if((!inNetwork.nodeExist(str[0]))||(str.length<=2))
//			{
//				continue;
//			}
//			else
//			{
//				Node source=inNetwork.getByName(str[0]);
//				for(int i=2;i<str.length;i++)
//				{
//					int k=i-1;
//					Node target_class=new Node("any",mn.Num_sample1,mn.Num_sample2);
//					target_class.deepcopyclass(k, source, inNetwork);
//					
//				}
//			}
//
//		}
//	    mappingreader.close();
	    ////////////////////////////////////////////////////////////////////////////
		
		//get confidence list
		ArrayList<Node> confidenceSet=new ArrayList<Node>();
		
		if(!mn.index_ifconf.equals("no"))
		{
			confidenceSet = fio.readConfidVectorfromFile(mn.inputfolder+mn.index_ifconf,inNetwork);
			for(Node i:confidenceSet)
			{
				i.setIsdif(true);
				for(Node inst : i.getContain())
				{
					inst.setIsdif(true);
				}
			}
		}
		//get start gene list
		
		ArrayList<Node> startPoint = fio.readStartVectorfromFile(mn.inputfolder+mn.start_gene,inNetwork);
		
		ArrayList<Node> endPoint=new ArrayList<Node>();
		
		//dir data_type sam1 sam2 keggnet ifconf start exp output topnum ifend ifppi
		//data/human_non_smoking_lung/ 2c_self_tfun1 60 60 output_net.txt no start.txt GSE19804_kegg_exp result50.txt 50 no no
		//data/mouse/ 2c_self_tfun1 9 9 output_net.txt no start.txt new_exp.txt 05252018_com_k 50 no no
		//data/mouse/ 2c_self_tfun1 9 9 output_net.txt no start.txt new_exp.txt 05252018_com_kp 50 no yes
		//data/mouse/ 2c_self_tfun1 3 3 output_net.txt no start.txt exp_control.txt 06152018_control_k 50 no no
		//data/mouse/ 2c_self_tfun1 3 3 output_net.txt no start_Casp4_fmod_has1_mir196b_col1a1_lum.txt exp_control.txt 03172019/control_kp 50 no ppi
		//data/yeast/cell_wall/ 2c_self_tfun1 6 6 output_net.txt no start.txt new_exp.txt result_newEXP //yeast_attr_map.txt, qvalue:1
		//data/yeast/cell_wall/ noexp_tfun3 6 6 output_net.txt no start.txt no result_dis 10 no ppi
		//data/paper_nsclc/ 2c_self_tfun1 60 60 output_net.txt no start.txt exp.txt result.txt
		//data/mouse/ 2c_self_tfun1 9 9 output_net.txt no start.txt new_exp.txt result.txt
		//data/lung_metastasis/ 2c_self_tfun1 659 23 output_net.txt no start_v2.txt exp.txt result.txt  //qvalue:1
		//data/ath/ 2c_self_tfun1 6 6 output_net.txt conf.txt start.txt fus3_exp.txt result_fus3  //qvalue:1, map_for_algo.txt
		//data/ChenYu_diabetes/ 2c_self_tfun1 11 10 output_net.txt no start_SLC5A2.txt exp13760.txt result_exp13760_slc5a2_50.txt 50 no yes
		//data/maize/V4/ 2c_self_tfun1 3 3 output_net.txt no start.txt tip_bif_exp.txt result_tip_bif_30 30 no ppi
		
		Pathway P=new Pathway();
	    //data/proteomic/BRCA2/ inte_pro 18 850 output_net.txt no start_NvP.txt BRCA_mrna_NvP_exp.txt result_50_NvP.txt 50 no yes BRCA_pro_NvP_exp.txt  
		//data/proteomic/BRCA/ inte_pro 850 5 output_net.txt no start_PvM.txt BRCA_mrna_PvM_exp2.txt result_50_PvM.txt 50 no yes BRCA_pro_PvM_exp2.txt
		//data/proteomic/BRCA2/ inte_pro 18 5 output_net.txt no start_PvM.txt BRCA_mrna_NvM_exp.txt result_30_NvM.txt 30 no yes BRCA_pro_NvM_exp.txt
		
		//data/proteomic/WashU/ inte_pro 4 4 output_net.txt no start_SPRY2_string.txt pro_exp_kegg result_50_MP.txt 50 no yes pro_exp_string
		//data/proteomic/WashU/ inte_pro 4 4 output_net.txt no start_SPRY2_string.txt exp_sprouty2suffi_vs_del_4vs4.txt result_50_MP2.txt 50 no yes pro_exp_string
		if(mn.dataType.equals("inte_pro"))
		{
			endPoint=fio.edgeWeight4Inte(mn.inputfolder+mn.index_ifexp,mn.inputfolder+mn.pro_exp_location, inNetwork, mn.Num_sample1,mn.Num_sample2,mn.Num_target, P);
			for(Node i:endPoint)
			{
				i.setIsdif(true);
				for(Node inst : i.getContain())
				{
					inst.setIsdif(true);
				}
			}
			mn.pv=P.calcuPathvalue(P.node_list);
		}
		

		//dir data_type sam1 sam2 keggnet ifconf start exp output topnum ifend ifppi pro_exp_location metab_exp num1meta num2meta
		//data/maize_metab/ inte_metab 5 5 output_net.txt no start_WWvWD_gene1.txt exp_gene_WWvsWD_a_5vs5.txt result_50_a2.txt 50 no no no exp_metab_WWvsWD_a_3vs3.txt 3 3
		if(mn.dataType.equals("inte_metab"))
		{
			endPoint=fio.edgeWeight4metab(mn.inputfolder+mn.index_ifexp,mn.inputfolder+mn.metab_exp, inNetwork, mn.Num_sample1,mn.Num_sample2,mn.Num_target, P, mn.num_sample1_metab, mn.num_sample2_metab);
			for(Node i:endPoint)
			{
				i.setIsdif(true);
				for(Node inst : i.getContain())
				{
					inst.setIsdif(true);
				}
			}
			mn.pv=P.calcuPathvalue(P.node_list);
		}
		
		
		
		//data/maize_drought/ 2c_self_tfun1 3 6 output_net.txt no start_WW_AvBC_gene1.txt WW_AvBC_exp.txt result_WW_AvBC_gene1_30 30 no ppi
		//data/maize_drought/ 2c_self_tfun1 9 9 output_net.txt no start_WWvWD_gene1.txt WWvWD_exp.txt result_WWvWD_gene1_30 30 no ppi
		//data/proteomic/BRCA/ 2c_self_tfun1 18 850 output_net.txt no start_NvP.txt BRCA_mRNA_NvP_exp.txt result_50_NvP_origin.txt 50 no yes
		//data/trupti_grant_20200909/ 2c_self_tfun1 5 5 output_net.txt no start_Foxc1.txt exp.txt result 50 no ppi no
		//data/trupti_grant_20200909/ 2c_self_tfun1 5 5 output_net.txt no start_Emx2.txt exp.txt result 50 no ppi no
		//data/trupti_grant_20200909/ 2c_self_tfun1 5 5 output_net.txt no start_foxo3.txt exp.txt result 50 no ppi no
		//data/trupti_grant_20200909/ 2c_self_tfun1 5 5 output_net.txt no start_Tcf21.txt exp.txt result 50 no ppi no
		//data/neptune/ 2c_self_tfun1 8 182 output_net.txt no start_sh3bp2.txt exp_glom.txt result_glom 50 no no no   
		//glom 8 182 (89MCD,93FSGS); n1g_n2g 6 97(47MCD,50FSGS); n1t_n2t 20 122(55MCD,67FSGS); tl 10 224(110MCD,114FSGS)
		//data/saintluke/ 2c_self_tfun1 4 4 output_net.txt no start_ERBB2.txt exp_squamous_TCGA_vs_normal_kegg result_squamous_tcga_vs_normal 50 no no no
		
		//data/proteomic/WashU/ 2c_self_tfun1 4 4 output_net.txt no start_SPRY2_mrna.txt pro_exp_kegg result_50_M.txt 50 no no no
		//data/proteomic/WashU/ 2c_self_tfun1 4 4 output_net.txt no start_SPRY2.txt exp_sprouty2suffi_vs_del_4vs4.txt result_50_M.txt 50 no ppi no 0 0 0
		//data/tarak_mouse/ 2c_self_tfun1 5 5 output_net.txt no start_Jak2.txt ctrl_vs_il6_5vs5_exp result_ctrl_il6_50.txt 50 no no no
		
		//data/Ron_soybean_flower/ 2c_self_tfun1 3 3 output_net.txt no start_ABA1.txt exp_c_vs_D.txt result_C_vs_D 50 no no no no 1 1
		//data/Ron_soybean_flower/ 2c_self_tfun1 3 3 output_net.txt conf.txt start_paper.txt exp_c_vs_D.txt result_C_vs_D_paper 50 no no no no 1 1
		
		//data/Hekmat_human/ 2c_self_tfun1 3 3 output_net.txt no start_cfl1.txt exp.txt result_cfl1_100 100 no no no no 1 1
		if((mn.dataType.equals("2c_self_tfun1"))&&(mn.index_ifend.equals("no")))
		{
			endPoint=fio.edgeWeight(mn.inputfolder+mn.index_ifexp, inNetwork, mn.Num_sample1,mn.Num_sample2,mn.Num_target, P);
			for(Node i:endPoint)
			{
				i.setIsdif(true);
				for(Node inst : i.getContain())
				{
					inst.setIsdif(true);
				}
			}
			mn.pv=P.calcuPathvalue(P.node_list);
		}
		else if((mn.dataType.equals("2c_self_tfun1"))&&(!mn.index_ifend.equals("no")))
		{
			fio.edgeWeight2(mn.inputfolder+mn.index_ifexp, inNetwork, mn.Num_sample1,mn.Num_sample2);
			endPoint=fio.readendVectorfromFile(mn.inputfolder+mn.index_ifend, inNetwork);
			P.getABCD(inNetwork, endPoint);
			mn.pv=P.calcuPathvalue(endPoint);
		}
		
		//data/2cancerCellline/ unisam_self_tfununi 1 1 output_net.txt no start_ptpn11.txt exp.txt result50 50 no no
		else if(mn.dataType.equals("unisam_self_tfununi"))
		{
			endPoint=fio.edgeWeight_NOreplicate(mn.inputfolder+mn.index_ifexp, inNetwork, mn.Num_target,P);
//			for(Node i: confidenceSet)
//			{
//				if(!endPoint.contains(i))
//				{
//					endPoint.add(i);
//				}
//			}
			for(Node i:endPoint)
			{
				i.setIsdif(true);
				for(Node inst : i.getContain())
				{
					inst.setIsdif(true);
				}
			}
			mn.pv=P.calcuPathvalue(P.node_list);
		}
		
		//data/mouse/ 2c_cuffdiff_tfun1 1 1 output_net.txt no start.txt il_kegg_diff.txt result_cuff_il.txt  //map_for_algo.txt
		//data/maize/ 2c_cuffdiff_tfun1 1 1 output_net.txt no start.txt bif2vsB73_base_kegg.txt result_cuff_bif2_base.txt
		else if(mn.dataType.equals("2c_cuffdiff_tfun1"))
		{
			endPoint=fio.edgeWeight_CuffDiff(mn.inputfolder+mn.index_ifexp, inNetwork,1.0);
			for(Node i:endPoint)
			{
				i.setIsdif(true);
				for(Node inst : i.getContain())
				{
					inst.setIsdif(true);
				}
			}
		}
	
		//dir data_type sam1 sam2 keggnet ifconf start exp output topnum ifend ifppi pro_exp_location metab_exp num1meta num2meta
		//data/Duolin_mouse/ dfgiven_tfun1 1 1 output_net.txt no start_Epha4.txt impress_log2fc_control_PI28 result_top50 50 no ppi
		//data/yeast/gal/ dfgiven_tfun1 1 1 output_net.txt no start.txt exp_gal80+gal.txt result_gal80+gal 50 no no
		//data/yeast/gal/ dfgiven_tfun1 1 1 output_net.txt conf.txt start_gal2_gal4_gal80.txt exp_gal80-gal.txt result_gal80-gal 20 no ppi
		//data/saintluke/ dfgiven_tfun1 4 4 output_net.txt no start_PIK3R2.txt DEG_squamous_cell_tcga_normal_vs_tumor.txt result_PIK3R2_tumor 50 no no no
		//data/saintluke/ dfgiven_tfun1 0 0 output_net.txt no start_ERBB2.txt DEG_given_tcga_normal_vs_own_minus_tumor_vs_own.txt result_200_ERBB2_given_tcga_normal_vs_own_minus_tumor_vs_own.txt 200 target_tcga_normal_vs_own_minus_tumor_vs_own_top200.txt no no no 0 0
		//data/saintluke/ dfgiven_tfun1 0 0 output_net.txt no start_ERBB2.txt DEG_given_tcga_normal_vs_tumor.txt result_200_ERBB2_given_tcga_normal_vs_tumor.txt 200 target_tcga_normal_vs_tumor_top200.txt no no no 0 0
		//data/saintluke/ dfgiven_tfun1 0 0 output_net.txt no start_ERBB2.txt DEG_given_tcga_tumor_vs_own.txt result_200_ERBB2_given_tcga_tumor_vs_own.txt 200 target_tcga_tumor_vs_own_top200.txt no no no 0 0
		//data/saintluke/ dfgiven_tfun1 0 0 output_net.txt no start_ERBB2.txt DEG_given_tcga_tumor_vs_own.txt result_ERBB2_given_tcga_tumor_vs_own.txt 0 target_tcga_tumor_vs_own.txt no no no 0 0
		//data/saintluke/crc/ dfgiven_tfun1 0 0 output_net.txt no start_spry2_kras.txt DEG_given_tcga_normal_vs_tumor.txt result_normal_vs_tumor 0 target_tcga_normal_vs_tumor_top200.txt no no no 0 0
		//data/saintluke/crc/ dfgiven_tfun1 0 0 output_net.txt no start_spry2_kras.txt DEG_given_tcga_normal_vs_own.txt result_normal_vs_own 0 target_tcga_normal_vs_own_top200.txt no no no 0 0
		else if((mn.dataType.equals("dfgiven_tfun1"))&&(mn.index_ifend.equals("no")))
		{
			endPoint=fio.edgeWeight_dfgiven(mn.inputfolder+mn.index_ifexp, inNetwork,mn.Num_target,P );
//			for(Node i: confidenceSet)
//			{
//				if(!endPoint.contains(i))
//				{
//					endPoint.add(i);
//				}
//			}
			for(Node i:endPoint)
			{
				i.setIsdif(true);
				for(Node inst : i.getContain())
				{
					inst.setIsdif(true);
				}
			}
			mn.pv=P.calcuPathvalue(P.node_list);
		}
		else if((mn.dataType.equals("dfgiven_tfun1"))&&(!mn.index_ifend.equals("no")))
		{
			fio.edgeWeight_dfgiven2(mn.inputfolder+mn.index_ifexp, inNetwork);
			endPoint=fio.readendVectorfromFile(mn.inputfolder+mn.index_ifend, inNetwork);
			P.getABCD(inNetwork, endPoint);
			mn.pv=P.calcuPathvalue(endPoint);
		}
		
		
		//data/yeast/cell_wall/ 2c_self_tfun1 6 6 output_net.txt no start.txt posrexp.txt result_rand 10 no ppi
		//data/yeast/cell_wall/ noexp_tfun3 6 6 output_net.txt no start.txt no result_dis 10 no ppi
		else if((mn.dataType.equals("noexp_tfun3"))&&(mn.index_ifend.equals("no")))
		{
			endPoint=fio.edgeWeight_noexp(inNetwork,mn.Num_target,P );
	//		endPoint=fio.edgeWeight_dfgiven(mn.inputfolder+mn.index_ifend, inNetwork, mn.Num_target, P);
			for(Node i:endPoint)
			{
				i.setIsdif(true);
				for(Node inst : i.getContain())
				{
					inst.setIsdif(true);
				}
			}
			mn.pv=P.calcuPathvalue(P.node_list);
		}
		else if((mn.dataType.equals("noexp_tfun3"))&&(!mn.index_ifend.equals("no")))
		{
			endPoint=fio.readendVectorfromFile(mn.inputfolder+mn.index_ifend, inNetwork);
			P.getABCD(inNetwork, endPoint);
			mn.pv=P.calcuPathvalue(endPoint);
		}
		
		
		
		else if(mn.dataType.equals("2c_list_tfun1"))
		{
			ArrayList<Node> JSlist=fio.edgeWeight(mn.inputfolder+mn.index_ifexp, inNetwork, mn.Num_sample1,mn.Num_sample2,mn.Num_target,P);
			endPoint=fio.readendVectorfromFile(mn.inputfolder+args[9], inNetwork);
		}
		//data/tarak_mouse/ time_list_tfun2 4 3 output_net.txt no start_Ctnnb1.txt exp.txt result_Ctnnb1_ppi 50 no ppi
		//run sample: data/yeast/osmotic/ time_list_tfun2 7 3 output_net.txt no start.txt exp.txt result.txt end.txt
		//data/yeast/osmotic/new_self/ time_list_tfun2 7 3 filtered_net.txt no start.txt exp.txt result01_noend_filter.txt 50 no no
		//data/yeast/osmotic/new_self/ time_list_tfun2 7 3 filtered_net.txt no start.txt exp.txt result01_filter.txt * end01.txt no 
		//data/yeast/osmotic/new_self/ time_list_tfun2 7 3 filtered_net.txt no start.txt exp.txt result01_noend30_filterPPI.txt 30 no ppi
		else if((mn.dataType.equals("time_list_tfun2"))&&(!mn.index_ifend.equals("no")))
		{
			
			fio.edgeWeight_time(mn.inputfolder+mn.index_ifexp, inNetwork, mn.Num_sample1*mn.Num_sample2);
			endPoint=fio.readendVectorfromFile(mn.inputfolder+args[10], inNetwork);
			P.getABCD(inNetwork, endPoint);
			mn.pv=P.calcuPathvalue(endPoint);
			
		}
		else if((mn.dataType.equals("time_list_tfun2"))&&(mn.index_ifend.equals("no")))
		{
			String userlistdir=fio.edgeWeight_time2(mn.inputfolder+mn.index_ifexp, inNetwork, mn.Num_sample1,mn.Num_sample2);
			
			endPoint=fio.readNendVectorfromFile(userlistdir+"/userlist", inNetwork,mn.Num_target);
			P.getABCD(inNetwork, endPoint);
			mn.pv=P.calcuPathvalue(endPoint);
		}
		else if(mn.dataType.equals("no_list_tfun2"))
		{
			endPoint=fio.readendVectorfromFile(mn.inputfolder+args[9], inNetwork);
		}
		else if(mn.dataType.equals("dream_time_tfunstep"))
		{
			inNetwork=new Network();
			inNetwork=fio.readNetworkfordream(mn.inputfolder+mn.netname, mn.Num_sample1,mn.Num_sample2);
			fio.edgeWeight_time(mn.inputfolder+mn.index_ifexp, inNetwork, mn.Num_sample1);
			endPoint=fio.readendVectorfromFile(mn.inputfolder+args[9], inNetwork);
			
			
		}
		//run sample: data/dream/bt20/ dream_tfundream 6 6 bignet.txt no start.txt 1.txt result1.txt end.txt
		else if(mn.dataType.equals("dream_tfundream"))
		{
			endPoint=fio.edgeWeight_dream(mn.inputfolder+mn.index_ifexp, inNetwork, mn.Num_sample1);
			for(Node i:endPoint)
			{
				i.setIsdif(true);
				for(Node inst : i.getContain())
				{
					inst.setIsdif(true);
				}
			}
		}
		
		
	//	P.getABCD(inNetwork, endPoint);
        
		
		
		

		
		
		ArrayList<Node> temEnd=new ArrayList<Node>();
		temEnd.addAll(0, endPoint);
		ArrayList<Node> temStart=new ArrayList<Node>();
		temStart.addAll(0, startPoint);
		ArrayList<Node> unfound=new ArrayList<Node>();
		unfound.addAll(0, endPoint);
		FindingPath fp=new FindingPath();
		
		//=====================================================================================================
//		double cover=fp.dijkstra_findmin(inNetwork, confidenceSet, temStart, temEnd, 1,mn.pv,mn.alpha);
//		System.out.println("start genes are: ");
//		for(Node j: temStart)
//		{
//			System.out.println(j.getname());
//		}
//		
//		if(cover<1)
//		{
//			unfound=mn.updateUnfound(unfound, temEnd);
//			
//			while(unfound.size()!=0)
//			{
//				System.out.println("============================in loop========================================");
//				System.out.println("unfound size= "+unfound.size() );
//				temEnd=mn.updateTemEnd(endPoint, temStart);
//				ArrayList<Node> newstart=new ArrayList<Node>();
//				newstart.add(unfound.get(0));
//				cover=fp.dijkstra_findmin(inNetwork, confidenceSet, newstart, temEnd, 1,mn.pv,mn.alpha);
//				System.out.println("start gene num: "+temStart.size()+", "+"end gene num: "+temEnd.size()+", "+"new start gene: "+newstart.get(0).getname());
//				
//				temStart=mn.updateTemStart(temStart, temEnd,newstart);
//				unfound=mn.updateUnfound(unfound, temEnd);
//				System.out.println("after updating, start genes are: ");
//				for(Node j: temStart)
//				{
//					System.out.println(j.getname());
//				}
//			}	
//		}
//		else
//		{
//			
//		}
//		
//		System.out.println("============================final========================================");	
//		temEnd=new ArrayList<Node>();
//		temEnd.addAll(0, endPoint);
		
		//=================================================statistic==================================================
//		int num_cla=0;
//		int num_inst=0;
//		ArrayList<Node> cla=new ArrayList<Node>();
//		ArrayList<Node> inst=new ArrayList<Node>();
//		for(Map.Entry<String, Node> i : inNetwork.getNodes().entrySet())
//		{
//			Node nd=i.getValue();
//			if((nd.getType().equals("abst"))&&(!(cla.contains(nd))))
//			{
//				
//				num_cla++;
//				cla.add(nd);
//				num_inst=num_inst+nd.getContain().size();
//			}
//			
//	
//		}
//		System.out.println("gene class num: "+num_cla);
//		System.out.println("gene inst num: "+num_inst);
//		
//		
//		
//		
//		File file1=new File(mn.inputfolder+"gold_internal.txt");
//		BufferedReader readerGI = new BufferedReader(new FileReader(file1));
//		File file2=new File(mn.inputfolder+"gold_tf.txt");
//		BufferedReader readerGTF = new BufferedReader(new FileReader(file2));
//		File file3=new File(mn.inputfolder+"map.txt");
//		BufferedReader map = new BufferedReader(new FileReader(file3));
//		File file4=new File(mn.inputfolder+"TF.txt");
//		BufferedReader readerTF = new BufferedReader(new FileReader(file4));
//		
//		String any="";
//		num_cla=0;
//		num_inst=0;
//		while ((any = readerGI.readLine()) != null)
//		{
//			if(inNetwork.nodeExist("abst#sce:"+any))
//			{
//				num_cla++;
//				Node nd=inNetwork.getByName("abst#sce:"+any);
//				num_inst=num_inst+nd.getContain().size();
//			}
//			else
//			{
//				System.out.println("gold internal gene missing: "+any);
//			}
//
//		}
//		System.out.println("gold internal gene class num: "+num_cla);
//		System.out.println("gold internal gene inst num: "+num_inst);
//		readerGI.close();
//		
//		any="";
//		num_cla=0;
//		num_inst=0;
//		while ((any = readerGTF.readLine()) != null)
//		{
//			if(inNetwork.nodeExist("abst#sce:"+any))
//			{
//				num_cla++;
//				Node nd=inNetwork.getByName("abst#sce:"+any);
//				num_inst=num_inst+nd.getContain().size();
//			}
//
//		}
//		System.out.println("gold TF gene class num: "+num_cla);
//		System.out.println("gold TF gene inst num: "+num_inst);
//		readerGTF.close();
//		
//		any="";
//		Hashtable<String, String> id2kegg = new Hashtable<String, String>();
//		while ((any = map.readLine()) != null)
//		{
//			String[] str=any.split("\t");
//			if(str.length==2)
//			{
//				id2kegg.put(str[1], str[0]);
//			}
//		}
//		
//		any="";
//		num_cla=0;
//		num_inst=0;
//		while ((any = readerTF.readLine()) != null)
//		{
//			
//			String name="abst#sce:"+id2kegg.get(any);
//			if(inNetwork.nodeExist(name))
//			{
//				num_cla++;
//				Node nd=inNetwork.getByName(name);
//				num_inst=num_inst+nd.getContain().size();
//			}
//			
//		}
//		System.out.println("TF gene class num: "+num_cla);
//		System.out.println("TF gene inst num: "+num_inst);
//		readerGTF.close();
//		
//		int atf=inNetwork.getByName("abst#sce:YMR172W").getContain().size()+inNetwork.getByName("abst#sce:YMR037C").getContain().size()+inNetwork.getByName("abst#sce:YKL062W").getContain().size()+inNetwork.getByName("abst#sce:YNL167C").getContain().size()-4;
//		int ainter=inNetwork.getByName("abst#sce:YLR113W").getContain().size()+inNetwork.getByName("abst#sce:YJL128C").getContain().size()+inNetwork.getByName("abst#sce:YLR362W").getContain().size()+inNetwork.getByName("abst#sce:YML004C").getContain().size()-4;
//		System.out.println("TF minus num: "+atf);
//		System.out.println("internal minus num: "+ainter);
		
//		ArrayList<String> pw=new ArrayList();
//		for(Node i : temEnd)
//		{
//			for(Node inst: i.getContain())
//			{
//				String p=inst.getPathway();
//				if(pw.contains(p))
//				{
//					continue;
//				}
//				else
//				{
//					pw.add(p);
//					System.out.println(p+" "+mn.pv.get(p));
//				}
//				
//			}
//		}
		
		
		
//		
//		for(Node i:temEnd)
//		{
//			System.out.println("Main: "+i.getname());
//		}
		//=============================================regular========================================================================
		if(mn.dataType.equals("2c_self_tfun1"))
		{
			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,1);
		}
		if(mn.dataType.equals("unisam_self_tfununi"))
		{
			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,1);
		}
		if(mn.dataType.equals("dfgiven_tfun1"))
		{
			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,1);
		}
		else if(mn.dataType.equals("time_list_tfun2"))
		{
			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,2);
		}
		else if(mn.dataType.equals("noexp_tfun3"))
		{
			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,3);
		}
		else if(mn.dataType.equals("inte_pro"))
		{
			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,1);
		}
		else if(mn.dataType.equals("inte_metab"))
		{
			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,1);
		}
		
		
	
		String temw = mn.inputfolder+mn.outname;
		Presentation pres=new Presentation();
//		pres.present_dream(temw, inNetwork);
		P=new Pathway();
//		pres.present(temw, inNetwork,P);
		pres.present_json2(temw, inNetwork, P, mn.inputfolder+"map_for_algo.txt",mn.inputfolder+"pathmapping.txt");
//		pres.present_json(temw, inNetwork, startPoint, confidenceSet, endPoint);
//		temw.close();
		
		Attribute atr=new Attribute();
		FileWriter atr1 = new FileWriter(mn.inputfolder+mn.outname+".attri");
		atr.gen_atr(inNetwork, atr1, temStart, endPoint,confidenceSet, mn.inputfolder+"map_for_algo.txt",mn.pv);
//		atr.gen_atr_kegg(inNetwork, atr1, temStart, endPoint, confidenceSet);
		atr1.close();
//		atr.tem_show(endPoint);
		
		///////////////////////////////////////////dream prior//////////////////////////////////////////////////////
//		Presentation pres=new Presentation();
////		FileWriter temw1 = new FileWriter(mn.inputfolder+"dismatrix.txt");
//		
//		pres.getMapping(mn.inputfolder+"map_bt20.txt");
//		FileWriter temw1 = new FileWriter(mn.inputfolder+"bt20sif.txt");
//		FileWriter temw2 = new FileWriter(mn.inputfolder+"bt20eda.txt");
//		for(Node e: endPoint)
//		{
//			temStart=new ArrayList<Node>();
//			temStart.add(e);
//			temEnd=new ArrayList<Node>();
//			temEnd.addAll(0, endPoint);
//			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,mn.alpha);
//
//			
////			pres.present_dream(temw1, endPoint);
//			pres.present_dream1(temw1, temw2, endPoint, e);
//			
//		}
//		temw1.close();
//		temw2.close();
		////////////////////////////////////////dream data////////////////////////////////////////////////////////////////
		
//		Presentation pres=new Presentation();
//		FileWriter temw1 = new FileWriter(mn.inputfolder+"subnet.txt");
//		pres.getMapping(mn.inputfolder+"map_bt20.txt");
//		for(Node e: endPoint)
//		{
//			System.out.println(e.getname());
//			temStart=new ArrayList<Node>();
//			temStart.add(e);
//			temEnd=new ArrayList<Node>();
//			temEnd.addAll(0, endPoint);
//			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,mn.alpha);
////			for(Node v: endPoint)
////			{
////				System.out.println(v.getValue());
////			}
//			pres.present_dream2(temw1, endPoint, e);
//		}
//		String name="abst#hsa:2260";
//		if(inNetwork.nodeExist(name))
//		{
//			Node e=inNetwork.getByName(name);
//			temStart=new ArrayList<Node>();
//			temStart.add(e);
//			temEnd=new ArrayList<Node>();
//			temEnd.addAll(0, endPoint);
//			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,mn.alpha);
//			pres.present_dream2(temw1, endPoint, e);
//		}
//		name="abst#hsa:2261";
//		if(inNetwork.nodeExist(name))
//		{
//			Node e=inNetwork.getByName(name);
//			temStart=new ArrayList<Node>();
//			temStart.add(e);
//			temEnd=new ArrayList<Node>();
//			temEnd.addAll(0, endPoint);
//			fp.dijkstra_findmin(inNetwork, confidenceSet, temStart,temEnd , 1,mn.pv,mn.alpha);
//			pres.present_dream2(temw1, endPoint, e);
//		}
//		
//		temw1.close();
		
		
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	}
	public ArrayList<Node> updateUnfound(ArrayList<Node> old, ArrayList<Node> newun)
	{
		ArrayList result=new ArrayList();
		for(Node i: newun)
		{
			if(old.contains(i))
			{
				result.add(i);
			}
		}
		return result;
	}
	public ArrayList<Node> updateTemEnd(ArrayList<Node> old, ArrayList<Node> start)
	{
		ArrayList<Node> result=new ArrayList<Node>();
		result.addAll(0, old);
		for(Node i: start)
		{
			if(!old.contains(i))
			{
				result.add(i);
			}
		}
		return result;
	}
	public ArrayList<Node> updateTemStart(ArrayList<Node> start, ArrayList<Node> End,ArrayList<Node> newstart)
	{
		ArrayList<Node> tem=new ArrayList<Node>();
		tem.addAll(0, start);
		for(Node i:start)
		{
			if(!End.contains(i))
			{
				tem.remove(i);
			}
		}
		for(Node i:newstart)
		{
			tem.add(i);
		}
		return tem;
	}
	
	
	
}
