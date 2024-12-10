package Algorithm;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map;








import org.rosuda.REngine.*;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;

public class Pathway {
	public Hashtable<String, Integer> pnum=new Hashtable<String, Integer>();
	public Hashtable<String, Integer> dpnum=new Hashtable<String, Integer>();
	public int abstnum=0;
	public int difnum=0;
	public ArrayList<String> pathlist=new ArrayList<String>();
	public Hashtable<String, Double> result=new Hashtable<String, Double>();
	public ArrayList<Node> node_list=new ArrayList<Node>();
	
	public double trans(double value, double base)
	{
		double v=Math.log(value);
		double b=Math.log(base);
		
		return v/b;
	}
	
	public void getABCD(Network nw, ArrayList<Node> al)
	{
		this.difnum=al.size();
		for(Map.Entry<String, Node> i : nw.getNodes().entrySet())
		{
			Node nd=i.getValue();
			if(nd.getType().equals("abst"))
			{
				this.abstnum++;
			}
			else
			{
				String pathway=nd.getPathway();
				
				if(this.pnum.containsKey(pathway))
				{
					int num=this.pnum.get(pathway)+1;
					this.pnum.put(pathway, num);
				}
				else
				{
					this.pnum.put(pathway,1);
				}
				
			}
		}
		for(int i=0;i<al.size();i++)
		{
			for(Node inst : al.get(i).getContain())
			{
				String pathway=inst.getPathway();
			
				if(!this.pathlist.contains(pathway))
				{
					this.pathlist.add(pathway);
				}
				if(this.dpnum.containsKey(pathway))
				{
					int num=this.dpnum.get(pathway)+1;
					this.dpnum.put(pathway, num);
				}
				else
				{
					this.dpnum.put(pathway,1);
				}
			}
		}
	}
	
	
	public Hashtable<String, Double> calcuPathvalue(ArrayList<Node> dif) throws IOException, REngineException, REXPMismatchException
	{
		RConnection con = new RConnection();
		
		
		for(Node i : dif)
		{		
				for(Node inst : i.getContain())
				{
					String pathname=inst.getPathway();
					if(pathname.equals("unknown"))
					{
						continue;
					}
					if(this.result.containsKey(pathname))
					{}
					else
					{
						int a=this.dpnum.get(pathname);
						int N=this.abstnum;
						int b=this.pnum.get(pathname)-a;
						int c=this.difnum-a;
						int d=N-a-b-c;
						int[] k={a,b,c,d};
						con.assign("a", k);
						con.assign("b", "greater");
						REXP yy=con.eval("fisher.test(matrix(a,2,2,byrow=TRUE),alternative = b)$p.value");
						double temp=1/-trans(yy.asDouble(),2);
						this.result.put(pathname, temp);
					}
				}			
		}
		con.close();
		return this.result;
		
	}

}
