package Algorithm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.rosuda.REngine.REXP;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.REngineException;
import org.rosuda.REngine.Rserve.RConnection;


public class Presentation {
	public Hashtable<String, String> kegg2id = new Hashtable<String, String>();
	public Hashtable<String, String> id2kegg = new Hashtable<String, String>();
	public Hashtable<String, String> pathname2realname = new Hashtable<String, String>();
	public Hashtable<String, Double> path2value=new Hashtable<String,Double>();
	public double maxfd=0.0;
	public double minfd=0.0;
	public Hashtable<String, String> name2link=new Hashtable<String,String>();
	
	public void read_info(String namemapdir,String pathnamedir) throws IOException
	{
		//map should be: id->kegg
		File mappingfile=new File(namemapdir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline;
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("\t");
			if(str.length>1)
			{
			this.kegg2id.put(str[1], str[0]);
			}
		}
		
		mappingfile=new File(pathnamedir);
		mappingreader=new BufferedReader(new FileReader(mappingfile));
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("  ",2);
			if(str.length>1)
			{
			//System.out.println(str[0]);
			//System.out.println(str[1]);
			this.pathname2realname.put(str[0], str[1]);
			}
		}
		
	
	}
	
	public String checktype(Node nd,ArrayList<Node> start, ArrayList<Node> conf, ArrayList<Node> end)
	{
		String type="";
		if(start.contains(nd.getBelong()))
			type="S";
		else if(conf.contains(nd.getBelong()))
			type="C";
		else if(end.contains(nd.getBelong()))
			type="E";
		return type;
	}
	
	public void present_json(FileWriter temw, Network net,ArrayList<Node> start, ArrayList<Node> conf, ArrayList<Node> end) throws IOException
	{
		StringBuilder nodes=new StringBuilder();
		StringBuilder edges=new StringBuilder();
		Hashtable<Node,Integer> check=new Hashtable<Node,Integer>();
		
		ArrayList<Node> leaf=new ArrayList<Node>();
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if((nd.getChildren().size()==0)&&(!nd.isIsdif())&&(!(nd.getType()=="compound"))&&(!(nd.getPrevious()==null)))
			{
				leaf.add(nd);
			}
		}
		while(!leaf.isEmpty())
		{
			Node del=leaf.get(0);
			Node del_p=del.getPrevious();
			del_p.getChildren().remove(del);
			if((del_p.getChildren().size()==0)&&(!del_p.isIsdif())&&(!(del_p.getType()=="compound"))&&(!(del_p.getPrevious()==null)))
			{
				leaf.add(del_p);
			}
			del.setPrevious(null);
			leaf.remove(del);
			
		}
		
		
		
		nodes.append("{\n\"nodes\":[\n");
		edges.append("\"edges\":[\n");
		int index=0;
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			
			if((!(nd.getPrevious()==null))&&(!(nd.getType().equals("abst"))))
			{
				String endnames=nd.getPrevious().getname()+nd.getname();
				Interaction interab=net.getByendnames(endnames);
				
				if(interab.getType().equals("type"))
				{
					if(nd.getPrevious().getPrevious()==null)
					{
						//
					}
					else
					{
						if((nd.getPrevious().getPrevious().getBelong()==nd.getBelong())&&(nd.getChildren().size()==0))
						{
							
						}
						else
						{
	
							if(!check.containsKey(nd))
							{
								String type=this.checktype(nd, start, conf, end);
								String pathway=nd.getPathway();
								check.put(nd, index++);
								String node="";
								if(index==1)
								{
									node="{ \"name\":\""+nd.getname()+"\",\"type\":\""+type+"\",\"pathway\":\""+pathway+"\"}";
								}
								else
								{
									node=",\n{ \"name\":\""+nd.getname()+"\",\"type\":\""+type+"\",\"pathway\":\""+pathway+"\"}";
								}
								
								nodes.append(node);
							}
							
							Node nd2=nd.getPrevious().getPrevious();
							if(!check.containsKey(nd2))
							{
								String type=this.checktype(nd2, start, conf, end);
								String pathway=nd2.getPathway();
								check.put(nd2, index++);
								String node=",\n{ \"name\":\""+nd2.getname()+"\",\"type\":\""+type+"\",\"pathway\":\""+pathway+"\"}";
								nodes.append(node);
							}
							
							
							if(index==2)
							{
								String edge="{\"source\":"+check.get(nd2)+",\"target\":"+check.get(nd)+",\"relation\":\"same gene\"}";
								edges.append(edge);
							}
							else
							{
								String edge=",\n{\"source\":"+check.get(nd2)+",\"target\":"+check.get(nd)+",\"relation\":\"same gene\"}";
								edges.append(edge);
							}
							
						}
						
					}
					
				}
				else
				{
					if(!check.containsKey(nd))
					{
						String type=this.checktype(nd, start, conf, end);
						String pathway=nd.getPathway();
						check.put(nd, index++);
						String node="";
						if(index==1)
						{
							node="{ \"name\":\""+nd.getname()+"\",\"type\":\""+type+"\",\"pathway\":\""+pathway+"\"}";
						}
						else
						{
							node=",\n{ \"name\":\""+nd.getname()+"\",\"type\":\""+type+"\",\"pathway\":\""+pathway+"\"}";
						}
						
						nodes.append(node);
					}
					Node nd2=nd.getPrevious();
					if(!check.containsKey(nd2))
					{
						String type=this.checktype(nd2, start, conf, end);
						String pathway=nd2.getPathway();
						check.put(nd2, index++);
						String node=",\n{ \"name\":\""+nd2.getname()+"\",\"type\":\""+type+"\",\"pathway\":\""+pathway+"\"}";
						nodes.append(node);
					}
					if(index==2)
					{
						String edge="{\"source\":"+check.get(nd2)+",\"target\":"+check.get(nd)+",\"relation\":\""+interab.getType()+"\"}";
						edges.append(edge);
					}
					else
					{
						String edge=",\n{\"source\":"+check.get(nd2)+",\"target\":"+check.get(nd)+",\"relation\":\""+interab.getType()+"\"}";
						edges.append(edge);
					}
					
				
				}
//			
				
				
			}
		}
		nodes.append("\n],\n");
		edges.append("\n]\n}");
		temw.write(nodes.toString()+edges.toString());
		
	temw.close();
		
	}
	
	public Node present(String temw, Network net, Pathway P) throws IOException, REngineException, REXPMismatchException
	{
		FileWriter fw=new FileWriter(temw+".tab");
		fw.write("ID1\tINTERACTION TYPE\tID2\n");
		FileWriter fw2=new FileWriter(temw+"_collapse.tab");
		fw2.write("gene name1\tINTERACTION TYPE\tgene name2\n");
		FileWriter ws_list=new FileWriter(temw+"_wslist.tab");
		Node virtual_node=new Node();
		virtual_node.setname("VirtualRoot");
		
		ArrayList<Node> leaf=new ArrayList<Node>();
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if((nd.getChildren().size()==0)&&(!nd.isIsdif())&&(!(nd.getType()=="compound"))&&(!(nd.getPrevious()==null)))
			{
				leaf.add(nd);
			}
		}
		while(!leaf.isEmpty())
		{
			Node del=leaf.get(0);
			Node del_p=del.getPrevious();
			del_p.getChildren().remove(del);
			if((del_p.getChildren().size()==0)&&(!del_p.isIsdif())&&(!(del_p.getType()=="compound"))&&(!(del_p.getPrevious()==null)))
			{
				leaf.add(del_p);
			}
			del.setPrevious(null);
			leaf.remove(del);
			
		}
		
		ArrayList<Node> al=new ArrayList<Node>();
		ArrayList<Node> al_inst=new ArrayList<Node>();
		Hashtable<String,Integer> pv=new Hashtable<String,Integer>();
		
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			
			if((!(nd.getPrevious()==null))&&(!(nd.getType().equals("abst"))))
			{
				
				String endnames=nd.getPrevious().getname()+nd.getname();
				Interaction interab=net.getByendnames(endnames);
				
				if(interab.getType().equals("type"))
				{
					if(nd.getPrevious().getPrevious()==null)
					{
						//
						al.add(nd.getBelong());
						al_inst.add(nd);
						if(pv.containsKey(nd.getPathway()))
						{
							pv.put(nd.getPathway(), pv.get(nd.getPathway())+1);
						}
						else
						{
							pv.put(nd.getPathway(), 1);
						}
						
						if(!nd.getChildren().isEmpty())
						{
							virtual_node.addTrueChildren(nd);
							nd.setTruefather(virtual_node);
							if(nd.getFc()>this.maxfd)
								this.maxfd=nd.getFc();
							if(nd.getFc()<this.minfd)
								this.minfd=nd.getFc();
						}
						
			//			System.out.println(nd.getname()+" "+nd.getFc()+" "+this.maxfd+" "+this.minfd);
					}
					else
					{
						if((nd.getPrevious().getPrevious().getBelong()==nd.getBelong())&&(nd.getChildren().size()==0))
						{
							
						}
						else
						{
							al.add(nd.getBelong());
							al_inst.add(nd);
							if(pv.containsKey(nd.getPathway()))
							{
								pv.put(nd.getPathway(), pv.get(nd.getPathway())+1);
							}
							else
							{
								pv.put(nd.getPathway(), 1);
							}
							fw.write(nd.getPrevious().getPrevious().getname()+"\tsame gene\t"+nd.getname());
							fw.write("\r\n");
							
			//				nd.getPrevious().getPrevious().addTrueChildren(nd);
			//				nd.setTruefather(nd.getPrevious().getPrevious());
							if(nd.getFc()>this.maxfd)
								this.maxfd=nd.getFc();
							if(nd.getFc()<this.minfd)
								this.minfd=nd.getFc();
							if(this.name2link.containsKey(nd.getname()))	
								this.name2link.put(nd.getname(), this.name2link.get(nd.getname())+",cross talk");
							else
								this.name2link.put(nd.getname(), "cross talk");
							if(this.name2link.containsKey(nd.getPrevious().getPrevious().getname()))
								this.name2link.put(nd.getPrevious().getPrevious().getname(), this.name2link.get(nd.getPrevious().getPrevious().getname())+",cross talk");
							else
								this.name2link.put(nd.getPrevious().getPrevious().getname(), "cross talk");
			//				System.out.println(nd.getname()+" "+nd.getFc()+" "+this.maxfd+" "+this.minfd);
						}
						
					}
					
				}
				//unmark from here
//				else if(interab.getType().equals("compound"))
//				{
//					al.add(nd.getBelong());
//					al_inst.add(nd);
//					if(pv.containsKey(nd.getPathway()))
//					{
//						pv.put(nd.getPathway(), pv.get(nd.getPathway())+1);
//					}
//					else
//					{
//						pv.put(nd.getPathway(), 1);
//					}
//					Node nd3=interab.getCompound();
//					fw.write(nd.getPrevious().getname()+"\t"+interab.getType()+"\t"+nd3.getname());
//					fw.write("\r\n");
//					fw.write(nd3.getname()+"\t"+interab.getType()+"\t"+nd.getname());
//					fw.write("\r\n");
//					String id1="";
//					String id2="";
//					id1=this.kegg2id.containsKey(nd.getPrevious().getname().split("#")[3])?this.kegg2id.get(nd.getPrevious().getname().split("#")[3]):nd.getPrevious().getname().split("#")[3];
//					id2=this.kegg2id.containsKey(nd.getname().split("#")[3])?this.kegg2id.get(nd.getname().split("#")[3]):nd.getname().split("#")[3];
//					fw2.write(id1+"\t"+interab.getType()+"\t"+id2+"\r\n");
//					
//					
//					nd.getPrevious().addTrueChildren(nd);
//					nd.setTruefather(nd.getPrevious());
//					if(nd.getFc()>this.maxfd)
//						this.maxfd=nd.getFc();
//					if(nd.getFc()<this.minfd)
//						this.minfd=nd.getFc();
//					for(int j=0;j<interab.getType().split("#").length;j++)
//					{
//						if(this.name2link.containsKey(nd.getname()))	
//							this.name2link.put(nd.getname(), this.name2link.get(nd.getname())+","+interab.getType().split("#")[j]);
//						else
//							this.name2link.put(nd.getname(), interab.getType().split("#")[j]);
//						if(this.name2link.containsKey(nd.getPrevious().getname()))
//							this.name2link.put(nd.getPrevious().getname(), this.name2link.get(nd.getPrevious().getname())+","+interab.getType().split("#")[j]);
//						else
//							this.name2link.put(nd.getPrevious().getname(), interab.getType().split("#")[j]);
//					}
//				}
				//unmark ends here
				else
				{
//					
					al.add(nd.getBelong());
					al_inst.add(nd);
					if(pv.containsKey(nd.getPathway()))
					{
						pv.put(nd.getPathway(), pv.get(nd.getPathway())+1);
					}
					else
					{
						pv.put(nd.getPathway(), 1);
					}
					fw.write(nd.getPrevious().getname()+"\t"+interab.getType()+"\t"+nd.getname());
					fw.write("\r\n");
					String id1="";
					String id2="";
					id1=this.kegg2id.containsKey(nd.getPrevious().getname().split("#")[3])?this.kegg2id.get(nd.getPrevious().getname().split("#")[3]):nd.getPrevious().getname().split("#")[3];
					id2=this.kegg2id.containsKey(nd.getname().split("#")[3])?this.kegg2id.get(nd.getname().split("#")[3]):nd.getname().split("#")[3];
					fw2.write(id1+"\t"+interab.getType()+"\t"+id2+"\r\n");
					
					
					nd.getPrevious().addTrueChildren(nd);
					nd.setTruefather(nd.getPrevious());
					if(nd.getFc()>this.maxfd)
						this.maxfd=nd.getFc();
					if(nd.getFc()<this.minfd)
						this.minfd=nd.getFc();
					for(int j=0;j<interab.getType().split("#").length;j++)
					{
						if(this.name2link.containsKey(nd.getname()))	
							this.name2link.put(nd.getname(), this.name2link.get(nd.getname())+","+interab.getType().split("#")[j]);
						else
							this.name2link.put(nd.getname(), interab.getType().split("#")[j]);
						if(this.name2link.containsKey(nd.getPrevious().getname()))
							this.name2link.put(nd.getPrevious().getname(), this.name2link.get(nd.getPrevious().getname())+","+interab.getType().split("#")[j]);
						else
							this.name2link.put(nd.getPrevious().getname(), interab.getType().split("#")[j]);
					}
					
			//		System.out.println(nd.getname()+" "+nd.getFc()+" "+this.maxfd+" "+this.minfd);
				}
//			
				
				
			}
		}
		for(Node nd:al_inst)
		{
			Node pnd=nd.getPrevious().getPrevious();
			if((pnd!=null)&&(pnd.getBelong()==nd.getBelong()))
			{
				pnd.addTruePathway(pnd.getPathway());
				pnd.addTruePathway(nd.getPathway());
				for(Node i:nd.getTrueChildren())
				{
					pnd.addTrueChildren(i);
					i.setTruefather(pnd);
				}
			}
		}
		
		P.getABCD(net, al);
		
//		System.out.println("Presentation.java: subnet size is: "+al.size());
		this.calcuPathvalue(al_inst, P, pv, this.path2value);
		
		//for ws_list
		ws_list.write("gene id\tgene name\tfold change\tpathway name\tpathway sig p-value\n");
		for(Node nd:al_inst)
		{
			String realname=this.kegg2id.containsKey(nd.getname().split("#")[3])?this.kegg2id.get(nd.getname().split("#")[3]):nd.getname().split("#")[3];
			System.out.print(nd.getPathway());
			String pname=this.pathname2realname.containsKey(nd.getPathway().substring(3))?this.pathname2realname.get(nd.getPathway().substring(3)):nd.getPathway();
			ws_list.write(nd.getname()+"\t"+realname+"\t"+nd.getFc()+"\t"+pname+"\t"+this.path2value.get(nd.getPathway()));
			ws_list.write("\r\n");
		}
		
	fw.close();
	fw2.close();
	ws_list.close();
	
	String[] flist={temw+".tab",temw+"_collapse.tab",temw+"_wslist.tab"};
	zip(flist,temw+".zip");
	return virtual_node;

	}
	
	public void zip(String[] filelist, String outputname) throws IOException {
        List<String> srcFiles = Arrays.asList(filelist);
        FileOutputStream fos = new FileOutputStream(outputname);
        ZipOutputStream zipOut = new ZipOutputStream(fos);
        for (String srcFile : srcFiles) {
            File fileToZip = new File(srcFile);
            FileInputStream fis = new FileInputStream(fileToZip);
            ZipEntry zipEntry = new ZipEntry(fileToZip.getName());
            zipOut.putNextEntry(zipEntry);
 
            byte[] bytes = new byte[1024];
            int length;
            while((length = fis.read(bytes)) >= 0) {
                zipOut.write(bytes, 0, length);
            }
            fis.close();
        }
        zipOut.close();
        fos.close();
    }
	
	public void dosomething(Node nd, FileWriter fw) throws IOException
	{
		double scalfc=0.0;
		if((nd.getFc()>0)&&(this.maxfd!=0.0))
			scalfc=nd.getFc()/this.maxfd;
		else if((nd.getFc()<0)&&(this.minfd!=0.0))
			scalfc=-nd.getFc()/this.minfd;
		
		
		String id="";
		id=this.kegg2id.containsKey(nd.getname().split("#")[3])?this.kegg2id.get(nd.getname().split("#")[3]):nd.getname().split("#")[3];
		fw.write("{\"name\":\""+id+"\",\r\n");
		
		fw.write("\"url\":\""+"http://www.genome.jp/dbget-bin/www_bget?"+nd.getname().split("#")[3]+"\",\r\n");
		
		fw.write("\"technos\":[");
//		System.out.println(nd.getname());
//		System.out.println(this.name2link.get(nd.getname()));
		int k=this.name2link.get(nd.getname()).split(",").length;
		for(int i=0;i<k-1;i++)
		{
			fw.write("\""+this.name2link.get(nd.getname()).split(",")[i]+"\",");
		}
		fw.write("\""+this.name2link.get(nd.getname()).split(",")[k-1]+"\"],\r\n");
		
		fw.write("\"satisfaction\":"+scalfc+",\r\n");
		
		k=nd.getTruePathway().size();
		if(k==2)
		{
			String p1=nd.getTruePathway().get(0);
			id=this.pathname2realname.containsKey(p1)?this.pathname2realname.get(p1):p1;
			fw.write("\"host\":{"+"\""+id+"\":[\"p-value:"+this.path2value.get(p1)+"\"],\r\n");
			String p2=nd.getTruePathway().get(1);
			id=this.pathname2realname.containsKey(p2)?this.pathname2realname.get(p2):p2;
			fw.write("\""+id+"\":[\"p-value:"+this.path2value.get(p2)+"\"]},\r\n");
		}
		else
		{
			String p1=nd.getPathway();
			id=this.pathname2realname.containsKey(p1)?this.pathname2realname.get(p1):p1;
			fw.write("\"host\":{"+"\""+id+"\":[\"p-value:"+this.path2value.get(p1)+"\"]},\r\n");
		}
			
		
		fw.write("\"children\":[");
		k=nd.getTrueChildren().size();
		if(k>0)
		{
			for(int i=0;i<k-1;i++)
			{
				Node cnd=nd.getTrueChildren().get(i);
				dosomething(cnd,fw);
				fw.write(",\r\n");
			}
			
			Node lcnd=nd.getTrueChildren().get(k-1);
			dosomething(lcnd,fw);
		}
		fw.write("]}\r\n");
		
		
		
	}
	
	public void present_json2(String fw, Network net, Pathway P,String mapdir,String pathmapdir) throws IOException, REngineException, REXPMismatchException
	{
		read_info(mapdir,pathmapdir);
		Node vr=present(fw,net,P);
//		System.out.println(this.maxfd+" "+this.minfd);
		FileWriter fw2=new FileWriter(fw+".json");
		fw2.write("{\"name\":\"VirtualRoot\",\r\n");
		fw2.write("\"children\":[\r\n");
		
		int k=vr.getTrueChildren().size();
		if(k>0)
		{
			for(int i=0;i<k-1;i++)
			{
				dosomething(vr.getTrueChildren().get(i),fw2);
				fw2.write(",\r\n");
			}
			dosomething(vr.getTrueChildren().get(k-1),fw2);
		}
		
		fw2.write("]}\r\n");
		fw2.close();
	}
	
	public void calcuPathvalue(ArrayList<Node> dif,Pathway P,Hashtable<String,Integer> ht, Hashtable<String, Double> pv) throws IOException, REngineException, REXPMismatchException
	{
		
		RConnection con = new RConnection();
		
		
		for(Node i : dif)
		{		
				
					String pathname=i.getPathway();
					if(pathname.equals("unknown"))
					{
						continue;
					}
					else
					{
						int a=ht.get(pathname);
						int N=P.abstnum;
						int b=P.pnum.get(pathname)-a;
						int c=dif.size()-a;
						int d=N-a-b-c;
						int[] k={a,b,c,d};
						con.assign("a", k);
						con.assign("b", "greater");
						REXP yy=con.eval("fisher.test(matrix(a,2,2,byrow=TRUE),alternative = b)$p.value");
						System.out.println("Presentation.java: "+ pathname+" "+yy.asDouble());
						if(!pv.containsKey(pathname))
						{
							pv.put(pathname, yy.asDouble());
						}
					}
						
		}
		
		con.close();
	}
	
	public void present_dream(FileWriter temw, Network net) throws IOException
	{
		ArrayList<Node> leaf=new ArrayList<Node>();
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if((nd.getChildren().size()==0)&&(!nd.isIsdif())&&(!(nd.getType()=="compound"))&&(!(nd.getPrevious()==null)))
			{
				leaf.add(nd);
			}
		}
		while(!leaf.isEmpty())
		{
			Node del=leaf.get(0);
			Node del_p=del.getPrevious();
			del_p.getChildren().remove(del);
			if((del_p.getChildren().size()==0)&&(!del_p.isIsdif())&&(!(del_p.getType()=="compound"))&&(!(del_p.getPrevious()==null)))
			{
				leaf.add(del_p);
			}
			del.setPrevious(null);
			leaf.remove(del);
			
		}
		
		
		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
		{
			Node nd=i.getValue();
			
			if((!(nd.getPrevious()==null))&&(!(nd.getType().equals("abst")))&&(nd.isIsdif()))
			{
				String endnames=nd.getPrevious().getname()+nd.getname();
				Interaction interab=net.getByendnames(endnames);
				
				if(interab.getType().equals("type"))
				{
					if(nd.getPrevious().getPrevious()==null)
					{
						//
					}
					else
					{
						if((nd.getPrevious().getPrevious().getBelong()==nd.getBelong())&&(nd.getChildren().size()==0))
						{
							
						}
						else
						{
							
//							temw.write(nd.getPrevious().getPrevious().getname()+"\tsame gene\t"+nd.getname());
//							temw.write("\r\n");
							Node target=nd.getBelong();
							Node p=this.find_dream(nd, target,0);
							if(p==null)
							{}
							else
							{
			//					double v=nd.getValue()-p.getValue();
								double v=-(nd.getValue()-p.getValue()-p.getIndex()*70)/(p.getIndex()*70);
								temw.write(p.getBelong().getname()+"\t"+nd.getBelong().getname()+"\t"+v);
								temw.write("\r\n");
							}
							
						}
						
					}
					
				}
				else
				{
//					
					
					Node target=nd.getBelong();
					Node p=this.find_dream(nd, target,0);
					if(p==null)
					{}
					else
					{
			//			double v=nd.getValue()-p.getValue();
						double v=-(nd.getValue()-p.getValue()-p.getIndex()*70)/(p.getIndex()*70);
						temw.write(p.getBelong().getname()+"\t"+nd.getBelong().getname()+"\t"+v);
						temw.write("\r\n");
					}
					
				}
//			
				
				
			}
		}
		
	temw.close();

	}
	
	public Node find_dream(Node nd, Node target,int num)
	{
		num++;
//		System.out.println(nd==null);
//		System.out.println(nd.getname()+"\t"+nd.getPrevious().getname());
		Node pnd= nd.getPrevious();
		boolean con3=(pnd==null);
		if(con3)
		{
			return pnd;
		}
		boolean con1=!pnd.isIsdif();
		boolean con2=(pnd.getType().equals("abst")||(pnd.getBelong()==target));
		
		
		while((con1||con2))
		{
		
			pnd=this.find_dream(pnd, target,num);
			break;
		}
		if(!(pnd==null))
		{
			pnd.setIndex(num);
		}
		
		return pnd;
	}
	
	public void present_dream1(FileWriter sif, FileWriter eda, ArrayList<Node> end, Node e) throws IOException
	{
		String id1=kegg2id.get(e.getname());
		String[] iid1=id1.split("\t");
		for(int i=0;i<iid1.length;i++)
		{
			for(Node v : end)
			{
				if((v.getValue()==Double.MAX_VALUE)||(v.getValue()==0))
				{}
				else
				{
					double d=(9-v.getValue())/8;
					String id2=kegg2id.get(v.getname());
					String[] iid2=id2.split("\t");
					for(int j=0;j<iid2.length;j++)
					{
						sif.write(iid1[i]+"\t1\t"+iid2[j]+"\r\n");
						
						eda.write(iid1[i]+" "+"(1)"+" "+iid2[j]+" = "+d+"\r\n");
					}
					
				}
			}
		}
		
	}
	
	public void present_dream2(FileWriter sub, ArrayList<Node> end, Node e) throws IOException  /////making subnet
	{
		String id1=kegg2id.get(e.getname());
		int l1=id1.split("\t").length;
	//	System.out.println(id1+" "+l1);
		if(l1==0)
		{
			for(Node v: end)
			{
				if((v.getValue()==Double.MAX_VALUE)||(v.getValue()==0))
				{}
				else
				{
					String id2=kegg2id.get(v.getname());
					int l2=id2.split("\t").length;
//					System.out.println(id2+" "+l2);
					for(int j=0;j<l2;j++)
					{
//						System.out.println("aaaaaaaaaaaaaaaaa");
						String name2=v.getname()+"?"+j;
						sub.write(e.getname()+"\tto\t"+name2+"\t"+v.getValue()+"\r\n");
					}
				}
			}
		}
		else
		{
		for(int i=0;i<l1;i++)
		{
			String name1=e.getname()+"?"+i;
			for(Node v: end)
			{
				if((v.getValue()==Double.MAX_VALUE)||(v.getValue()==0))
				{}
				else
				{
					
					String id2=kegg2id.get(v.getname());
					int l2=id2.split("\t").length;
//					System.out.println(id2+" "+l2);
					for(int j=0;j<l2;j++)
					{
//						System.out.println("aaaaaaaaaaaaaaaaa");
						String name2=v.getname()+"?"+j;
						sub.write(name1+"\tto\t"+name2+"\t"+v.getValue()+"\r\n");
					}
				}
			}
		}
		}
	}
	
	public void getMapping(String dir) throws IOException
	{
		
		File mappingfile=new File(dir);
		BufferedReader mappingreader=new BufferedReader(new FileReader(mappingfile));
		String mappingline=null;
		while((mappingline=mappingreader.readLine())!=null)
		{
			String[] str=mappingline.split("\t",2);
			
			this.kegg2id.put(str[0], str[1]);
			String[] sub=str[1].split("\t");
			int l=sub.length;
			for(int i=0;i<l;i++)
			{
				this.id2kegg.put(sub[i], str[0]);
			}
				
		}
						
		mappingreader.close();
		
	}

//	public ArrayList<Path> sortIncrement(ArrayList<Path> al)
//	{
//		Collections.sort(al,new Comparator<Path>()
//		{
//			public int compare(Path arg0, Path arg1) 
//			{
//				double tem=(arg0.getValue() - arg1.getValue())*10000000;
//				return (int)tem;
//            }
//        });
//		return al;
//	}
//	public void showPath(ArrayList<Path> paths, String filename, Network net, double percent, String nettype, boolean ifper, int keepnum) throws IOException
//	{
//		paths=sortIncrement(paths);
//		int per=(int)(paths.size()*percent);
//		if(per<1)
//			per=1;
//		if(!ifper)
//			per=keepnum;
//		List<Path> subpaths=paths.subList(0, per);
//		
//		String output_net=filename+"net.txt";
//		String output_attri=filename+"attri.txt";
//		String output_info=filename+"info.txt";
//		FileWriter fw_net=new FileWriter(output_net);
//		FileWriter fw_attri=new FileWriter(output_attri);
//		FileWriter fw_info=new FileWriter(output_info);
//		
//		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
//		{
//			if((i.getValue().getType().equals("gene"))||(i.getValue().getType().equals("group")))
//			{
//			fw_attri.write(i.getValue().getname()+" "+i.getValue().getFlag()+" "+i.getValue().getPathway()+" "+i.getValue().getBelong().getname());
//			fw_attri.write("\r\n");
//			}
//		}
//		
//		for(int i=0;i<subpaths.size();i++)
//		{
//			Path current=subpaths.get(i);
//			int length=current.getNodes().size();
//			int line=0;
//			if(length==1)
//			{
//				System.out.println("Does this case exist?");
//			}
//			else 
//			{
//				ArrayList<Node> tem = new ArrayList<Node>();
//				for(int z=0;z<current.getNodes().size();z++)
//				{
//					tem.add(current.getNodes().get(z));
//				}
//				
//				tem.remove(length-1);
//				tem.remove(0);
//				int newlength=length-2;
//				for(int j=newlength-1;j>=1;j--)
//				{
//					Node back0=tem.get(j);
//					Node back1=tem.get(j-1);
//					String str01=back0.getname()+back1.getname();
//					Interaction inter01=net.getByendnames(str01);
//					if(inter01.getType().equals("type"))
//					{
//						Node back2=tem.get(j-2);
//						fw_net.write(back0.getname()+" <same class> "+back2.getname());
//						fw_net.write("\r\n");
//						line++;
//						j--;
//					}
//					else if(inter01.getType().equals("compound"))
//					{
//						if(nettype.equals("gc"))
//						{
//							fw_net.write(back0.getname()+" <compound> "+inter01.getCompound().getname());
//							fw_net.write("\r\n");
//							line++;
//							fw_net.write(inter01.getCompound().getname()+" <compound> "+back1.getname());
//							fw_net.write("\r\n");
//							line++;
//						}
//						if(nettype.equals("gg"))
//						{
//							fw_net.write(back0.getname()+" <compound> "+back1.getname());
//							fw_net.write("\r\n");
//							line++;
//						}
//						
//					}
//					else
//					{
//						fw_net.write(back0.getname()+" <"+inter01.getType()+"> "+back1.getname());
//						fw_net.write("\r\n");
//						line++;
//					}
//					
//				}
//				
//			}
//			fw_info.write(line+" "+current.getValue()+" "+current.getConfinum());
//			fw_info.write("\r\n");
//			
//		}
//		fw_attri.close();
//		fw_info.close();
//		fw_net.close();
//	}
//
//	public void showPath2(BufferedReader rd, int pathnum, String filename, Network net, double percent, String nettype, boolean ifper, int keepnum) throws IOException
//	{
//		
//		int per=(int)(pathnum*percent);
//		if(per<1)
//			per=1;
//		if(!ifper)
//			per=keepnum;
//		ArrayList<Path> subpaths=new ArrayList<Path>();
//		String temString=null;
//		while((temString=rd.readLine())!=null)
//		{
//			Path a=new Path();
//			a.setfromarray(temString, net);
//			if(subpaths.size()==per)
//			{
//				if(a.getValue()>=subpaths.get(per-1).getValue())
//				{
//					continue;
//				}
//				else 
//				{
//					subpaths.remove(per-1);
//					subpaths.add(a);
//					subpaths=sortIncrement(subpaths);
//				}
//			}
//			else 
//			{
//				subpaths.add(a);
//				if(subpaths.size()==per)
//				{
//					subpaths=sortIncrement(subpaths);
//				}
//			}
//		}
//
//		
//		String output_net=filename+"net.txt";
//		String output_attri=filename+"attri.txt";
//		String output_info=filename+"info.txt";
//		FileWriter fw_net=new FileWriter(output_net);
//		FileWriter fw_attri=new FileWriter(output_attri);
//		FileWriter fw_info=new FileWriter(output_info);
//		
//		for(Map.Entry<String, Node> i : net.getNodes().entrySet())
//		{
//			if((i.getValue().getType().equals("gene"))||(i.getValue().getType().equals("group")))
//			{
//			fw_attri.write(i.getValue().getname()+" "+i.getValue().getFlag()+" "+i.getValue().getPathway()+" "+i.getValue().getBelong().getname());
//			fw_attri.write("\r\n");
//			}
//		}
//		
//		for(int i=0;i<subpaths.size();i++)
//		{
//			Path current=subpaths.get(i);
//			int length=current.getNodes().size();
//			int line=0;
//			if(length==1)
//			{
//				System.out.println("Does this case exist?");
//			}
//			else 
//			{
//				ArrayList<Node> tem = new ArrayList<Node>();
//				for(int z=0;z<current.getNodes().size();z++)
//				{
//					tem.add(current.getNodes().get(z));
//				}
//				
//				tem.remove(length-1);
//				tem.remove(0);
//				int newlength=length-2;
//				for(int j=newlength-1;j>=1;j--)
//				{
//					Node back0=tem.get(j);
//					Node back1=tem.get(j-1);
//					String str01=back0.getname()+back1.getname();
//					Interaction inter01=net.getByendnames(str01);
//					if(inter01.getType().equals("type"))
//					{
//						Node back2=tem.get(j-2);
//						fw_net.write(back0.getname()+" <same class> "+back2.getname());
//						fw_net.write("\r\n");
//						line++;
//						j--;
//					}
//					else if(inter01.getType().equals("compound"))
//					{
//						if(nettype.equals("gc"))
//						{
//							fw_net.write(back0.getname()+" <compound> "+inter01.getCompound().getname());
//							fw_net.write("\r\n");
//							line++;
//							fw_net.write(inter01.getCompound().getname()+" <compound> "+back1.getname());
//							fw_net.write("\r\n");
//							line++;
//						}
//						if(nettype.equals("gg"))
//						{
//							fw_net.write(back0.getname()+" <compound> "+back1.getname());
//							fw_net.write("\r\n");
//							line++;
//						}
//						
//					}
//					else
//					{
//						fw_net.write(back0.getname()+" <"+inter01.getType()+"> "+back1.getname());
//						fw_net.write("\r\n");
//						line++;
//					}
//					
//				}
//				
//			}
//			fw_info.write(line+" "+current.getValue()+" "+current.getConfinum());
//			fw_info.write("\r\n");
//			
//		}
//		fw_attri.close();
//		fw_info.close();
//		fw_net.close();
//	}
//	public void showPath2Json(ArrayList<Path> paths, String file, Network net, double percent, String nettype, boolean ifper, int keepnum) throws IOException
//	{
//		paths=sortIncrement(paths);
//		int per=(int)(paths.size()*percent);
//		if(per<1)
//			per=1;
//		if(!ifper)
//			per=keepnum;
//		List<Path> subpaths=paths.subList(0, per);
//		
//		Hashtable<Node, Integer> n2i=new Hashtable<Node, Integer>();
//		Hashtable<Integer, Node> i2n=new Hashtable<Integer, Node>();
//		FileWriter writer=new FileWriter(file);
//		int nodeNumber=0;
////		StringBuilder pathString = new StringBuilder();
//		File temfile=new File(file+"tem");
//		FileWriter temwriter=new FileWriter(temfile);
//		
////		pathString.append("\n],\"edge\":[\n");
//		temwriter.write("\n],\"edge\":[\n");
//		
//		for(int i=0;i<subpaths.size();i++)
//		{
//		
//			Path current=subpaths.get(i);
//			int line=0;
//			int length=current.getNodes().size();
//			
//			if(length==1)
//			{
//				System.out.println("I dont't think this case exists!");
//			}
//			else 
//			{
//				ArrayList<Node> tem = new ArrayList<Node>();
//				for(int z=0;z<current.getNodes().size();z++)
//				{
//					tem.add(current.getNodes().get(z));
//				}
//				
//				tem.remove(length-1);
//				tem.remove(0);
//				int newlength=length-2;
//				for(int j=newlength-1;j>=1;j--)
//				{
//					Node back0=tem.get(j);
//					Node back1=tem.get(j-1);
//					String str01=back0.getname()+back1.getname();
//					Interaction inter01=net.getByendnames(str01);
//					
//					if(!n2i.containsKey(back0))
//					{
//						n2i.put(back0, nodeNumber);
//						i2n.put(nodeNumber, back0);
//						nodeNumber++;
//					}
//					if((!n2i.containsKey(back1))&&(!back1.getType().equals("abst")))
//					{
//						n2i.put(back1, nodeNumber);
//						i2n.put(nodeNumber, back1);
//						nodeNumber++;
//					}
//
//					
//					if(inter01.getType().equals("type"))
//					{
//						Node back2=tem.get(j-2);
//						if(!n2i.containsKey(back2))
//						{
//							n2i.put(back2, nodeNumber);
//							i2n.put(nodeNumber, back2);
//							nodeNumber++;
//						}
//		//				pathString.append("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back2)+",\"relation\":\"same_class\"},\n");
//						temwriter.write("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back2)+",\"relation\":\"same_class\"},\n");
//						
//						line++;
//						j--;
//					}
//					else if(inter01.getType().equals("compound"))
//					{
//						if(nettype.equals("gc"))
//						{
//							Node com=inter01.getCompound();
//							if(!n2i.containsKey(com))
//							{
//								n2i.put(com, nodeNumber);
//								i2n.put(nodeNumber, com);
//								nodeNumber++;
//							}
//			//				pathString.append("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(com)+",\"relation\":\""+inter01.getType()+"\"},\n");
//							temwriter.write("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(com)+",\"relation\":\""+inter01.getType()+"\"},\n");
//							
//							line++;
//			//				pathString.append("{\"source\":"+n2i.get(com)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"},\n");
//							temwriter.write("{\"source\":"+n2i.get(com)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"},\n");
//							
//							line++;
//						}
//						if(nettype.equals("gg"))
//						{
//							temwriter.write("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"},\n");
//							line++;
//						}
//						
//					}
//					else
//					{
//		//				pathString.append("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"},\n");
//						temwriter.write("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"},\n");
//						
//						line++;
//					}
//					
//				}
//				
//			}
//		}
////		pathString.append("]\n}\n");
//		temwriter.write("]\n}\n");
//		temwriter.close();
//		
//		StringBuilder nodeString=new StringBuilder();
//		nodeString.append("{\n\"nodes\":[\n");
//		for(int i=0;i<nodeNumber;i++)
//		{
//			nodeString.append("{ \"name\":\""+i2n.get(i).getname()+"\",\"type\":\""+i2n.get(i).getFlag()+"\",\"pathway\":\""+i2n.get(i).getPathway()+"\"},\n");
//		}
//		writer.write(nodeString.toString());
//		
////		File rfile=new File(file+"tem");
//		BufferedReader reader = new BufferedReader(new FileReader(temfile));
//		String tempString = null;
//    	while ((tempString = reader.readLine()) != null)
//    	{
//    		writer.write(tempString);
//    		writer.write("\r\n");
//    	}
//		
////		writer.write(pathString.toString());
//		writer.close();
//	//	File f=new File(file+"tem");
//		temfile.delete();
//	}
//	public void showPath2Json2(BufferedReader rd, int pathnum, String file, Network net, double percent, String nettype, boolean ifper, int keepnum) throws IOException
//	{		
//		int per=(int)(pathnum*percent);
//		if(per<1)
//			per=1;
//		if(!ifper)
//			per=keepnum;
//		ArrayList<Path> subpaths=new ArrayList<Path>();
//		String temString=null;
//		while((temString=rd.readLine())!=null)
//		{
//			Path a=new Path();
//			a.setfromarray(temString, net);
//			if(subpaths.size()==per)
//			{
//				if(a.getValue()>=subpaths.get(per-1).getValue())
//				{
//					continue;
//				}
//				else 
//				{
//					subpaths.remove(per-1);
//					subpaths.add(a);
//					subpaths=sortIncrement(subpaths);
//				}
//			}
//			else 
//			{
//				subpaths.add(a);
//				if(subpaths.size()==per)
//				{
//					subpaths=sortIncrement(subpaths);
//				}
//			}
//		}
//		
//		
//		
//		
//		Hashtable<Node, Integer> n2i=new Hashtable<Node, Integer>();
//		Hashtable<Integer, Node> i2n=new Hashtable<Integer, Node>();
//		FileWriter writer=new FileWriter(file);
//		int nodeNumber=0;
////		StringBuilder pathString = new StringBuilder();
//		File temfile=new File(file+"tem");
//		FileWriter temwriter=new FileWriter(temfile);
//		
////		pathString.append("\n],\"edge\":[\n");
//		temwriter.write("],\n\"edges\":[\n");
//		
//		for(int i=0;i<subpaths.size();i++)
//		{
//		
//			Path current=subpaths.get(i);
//			int line=0;
//			int length=current.getNodes().size();
//			
//			if(length==1)
//			{
//				System.out.println("I dont't think this case exists!");
//			}
//			else 
//			{
//				ArrayList<Node> tem = new ArrayList<Node>();
//				for(int z=0;z<current.getNodes().size();z++)
//				{
//					tem.add(current.getNodes().get(z));
//				}
//				
//				tem.remove(length-1);
//				tem.remove(0);
//				int newlength=length-2;
//				for(int j=newlength-1;j>=1;j--)
//				{
//					Node back0=tem.get(j);
//					Node back1=tem.get(j-1);
//					String str01=back0.getname()+back1.getname();
//					Interaction inter01=net.getByendnames(str01);
//					
//					if(!n2i.containsKey(back0))
//					{
//						n2i.put(back0, nodeNumber);
//						i2n.put(nodeNumber, back0);
//						nodeNumber++;
//					}
//					if((!n2i.containsKey(back1))&&(!back1.getType().equals("abst")))
//					{
//						n2i.put(back1, nodeNumber);
//						i2n.put(nodeNumber, back1);
//						nodeNumber++;
//					}
//
//					
//					if(inter01.getType().equals("type"))
//					{
//						Node back2=tem.get(j-2);
//						if(!n2i.containsKey(back2))
//						{
//							n2i.put(back2, nodeNumber);
//							i2n.put(nodeNumber, back2);
//							nodeNumber++;
//						}
//		//				pathString.append("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back2)+",\"relation\":\"same_class\"},\n");
//						temwriter.write("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back2)+",\"relation\":\"same_class\"}");
//						if((j-2==0)&&(i==subpaths.size()-1))
//							temwriter.write("\n");
//						else {
//							temwriter.write(",\n");
//						}
//						
//						line++;
//						j--;
//					}
//					else if(inter01.getType().equals("compound"))
//					{
//						if(nettype.equals("gc"))
//						{
//							Node com=inter01.getCompound();
//							if(!n2i.containsKey(com))
//							{
//								n2i.put(com, nodeNumber);
//								i2n.put(nodeNumber, com);
//								nodeNumber++;
//							}
//			//				pathString.append("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(com)+",\"relation\":\""+inter01.getType()+"\"},\n");
//							temwriter.write("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(com)+",\"relation\":\""+inter01.getType()+"\"},\n");
//							
//							line++;
//			//				pathString.append("{\"source\":"+n2i.get(com)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"},\n");
//							temwriter.write("{\"source\":"+n2i.get(com)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"}");
//							if((j-1==0)&&(i==subpaths.size()-1))
//								temwriter.write("\n");
//							else {
//								temwriter.write(",\n");
//							}
//							line++;
//						}
//						if(nettype.equals("gg"))
//						{
//							temwriter.write("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"}");
//							if((j-1==0)&&(i==subpaths.size()-1))
//								temwriter.write("\n");
//							else {
//								temwriter.write(",\n");
//							}
//							line++;
//						}
//						
//					}
//					else
//					{
//		//				pathString.append("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"},\n");
//						temwriter.write("{\"source\":"+n2i.get(back0)+",\"target\":"+n2i.get(back1)+",\"relation\":\""+inter01.getType()+"\"}");
//						if((j-1==0)&&(i==subpaths.size()-1))
//							temwriter.write("\n");
//						else {
//							temwriter.write(",\n");
//						}
//						line++;
//					}
//					
//				}
//				
//			}
//		}
////		pathString.append("]\n}\n");
//		temwriter.write("]\n}\n");
//		temwriter.close();
//		
//		StringBuilder nodeString=new StringBuilder();
//		nodeString.append("{\n\"nodes\":[\n");
//		for(int i=0;i<nodeNumber;i++)
//		{
//			nodeString.append("{ \"name\":\""+i2n.get(i).getname()+"\",\"type\":\""+i2n.get(i).getFlag()+"\",\"pathway\":\""+i2n.get(i).getPathway()+"\"}");
//			if(i!=nodeNumber-1)
//				nodeString.append(",\n");
//			else {
//				nodeString.append("\n");
//			}
//		}
//		writer.write(nodeString.toString());
//		
////		File rfile=new File(file+"tem");
//		BufferedReader reader = new BufferedReader(new FileReader(temfile));
//		String tempString = null;
//    	while ((tempString = reader.readLine()) != null)
//    	{
//    		writer.write(tempString);
//    		writer.write("\r\n");
//    	}
//		reader.close();
////		writer.write(pathString.toString());
//		writer.close();
//	//	File f=new File(file+"tem");
//		temfile.delete();
//	}
}
