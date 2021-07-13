import java.io.*;
import java.lang.*;
import java.util.*;


public class countDMScodons
{	
	public int[][] codoncount;
	public int[][] aacount;
	public String[] CODONS = new String[] {"AAA", "AAG", "AAC", "AAT",
											"AGA", "AGG", "AGC", "AGT", 
											"ACA", "ACG", "ACC","ACT", 
											"ATA", "ATG", "ATC","ATT", 
											"GAA", "GAG", "GAC","GAT", 
											"GGA", "GGG", "GGC","GGT", 
											"GCA", "GCG", "GCC","GCT",
											"GTA", "GTG", "GTC","GTT",
											"CAA", "CAG", "CAC","CAT",
											"CGA", "CGG", "CGC","CGT",
											"CCA", "CCG", "CCC","CCT",
											"CTA", "CTG", "CTC","CTT",
											"TAA", "TAG", "TAC","TAT",
											"TGA", "TGG", "TGC","TGT",
											"TCA", "TCG", "TCC","TCT",
											"TTA", "TTG", "TTC","TTT"};
											
	public String[] AACIDS = new String[] {"K", "N", "R", "S", "T", "I", "M", "E","D","G","A","V","Q","H","P","L","*","Y","W","C","F"};

	
	public static void main(String args[])
	{
		new countDMScodons(args[0], Integer.valueOf(args[1]), Integer.valueOf(args[2]));
	}
	
	public countDMScodons(String filename, int start, int stop)
	{
		 codoncount = new int[(stop-start)/3+1][64];
		 aacount = new int[(stop-start)/3+1][21];
		 
		 for(int i =0; i < codoncount.length;i++)
		 	 for(int j =0; j <codoncount[i].length; j++)
		 	 	codoncount[i][j] = 0;
		 	
		 for(int i =0; i < aacount.length;i++)
		 	 for(int j =0; j < aacount[i].length; j++)
		 	 	aacount[i][j] = 0;
		 
		 readSAM(filename, start, stop);
	}
	
	public void readSAM(String filename, int start, int stop)
	{
		try
		{
			BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
			
			String temp = br.readLine(); 
			
			while (temp != null)
			{
				if (temp.charAt(0) == '@')
				{
					temp = br.readLine(); //move through header
				}
				else
				{
					StringTokenizer tk1 = new StringTokenizer(temp, "\t");
					String name1 = tk1.nextToken(); //read name
					
					temp = br.readLine(); //read next record (should be mate)
					
					StringTokenizer tk2 = new StringTokenizer(temp, "\t");
					String name2 = "";
					if (tk2.hasMoreTokens())
						 name2 = tk2.nextToken(); //read name
					
					if(name2.equals(name1))
					{						
						SAMString one = new SAMString(name1, tk1.nextToken(), tk1.nextToken(), Integer.valueOf(tk1.nextToken()), Integer.valueOf(tk1.nextToken()), tk1.nextToken(), tk1.nextToken(), Integer.valueOf(tk1.nextToken()), Integer.valueOf(tk1.nextToken()), tk1.nextToken(), tk1.nextToken()); //name, flag, genome_name, start, MQ, cigar, mate,  of mate1
						SAMString two = new SAMString(name2, tk2.nextToken(), tk2.nextToken(), Integer.valueOf(tk2.nextToken()), Integer.valueOf(tk2.nextToken()), tk2.nextToken(), tk2.nextToken(), Integer.valueOf(tk2.nextToken()), Integer.valueOf(tk2.nextToken()), tk2.nextToken(), tk2.nextToken()); //name, flag, genome_name, start, MQ, cigar, mate,  of mate1
						
						countCodons(start, stop, new SAMJoin(one, two));
						
						
						temp = br.readLine(); //read next record
					}
					else 
					{
						System.out.println(name1 + " has no mate or SAM needs to be sorted");
					}
				}
			}
			br.close();
		}
		catch (Exception e)
		{
			System.out.println("crap");
			e.printStackTrace();
		}
		writeCodons(filename, start, stop);
		writeAA(filename, start, stop);
	}
	
	public void countCodons(int start, int stop, SAMJoin cluster)
	{
		//System.out.print(cluster.getName()+"\t");
		
		SAMString mate1 = cluster.getMate1();
		SAMString mate2 = cluster.getMate2();
		String full = "";
		
		if (start > cluster.getStart() && stop < cluster.getEnd())
		{
			for (int i = 0; (start-cluster.getStart()+i+3) < (stop-cluster.getStart()); i = i+3)
			{
				try
				{
					if ((start-cluster.getStart()+i+3) < mate1.getaRead().length() && i < (mate2.getStart()-mate1.getStart())) //if mate 1 does not overlap mate 2
					{
						full = full + mate1.getaRead().substring(start-cluster.getStart()+i, start-cluster.getStart()+i+3);
					}
					else if((start-cluster.getStart()+i+3) < mate1.getaRead().length() && i +start > mate2.getStart() && (start+i-mate2.getStart()+3) < mate2.getaRead().length()) //if still within mate1 and mate2 has started and hasn't ended (i.e. overlapped)
					{							
						String s = "";
						for (int j = i; j < i+3; j++)
						{
							int offset = start+j-mate2.getStart();
							if (mate1.getQual().charAt(start-cluster.getStart()+j) > mate2.getQual().charAt(offset)) //(check who has better quality at this base, add to small substring
								s = s + mate1.getaRead().charAt(start-cluster.getStart()+j);
							else
								s = s + mate2.getaRead().charAt(offset);
						}
						
						full = full + s; //add highest quality substring to  full
					}
					else if(i +start > mate2.getStart() && start+i-mate2.getStart()+3 < mate2.getaRead().length()) //mate 2 non-overlap
					{
						int offset = start+i-mate2.getStart();
						full = full + mate2.getaRead().substring(offset,offset+3);
					}	
					
				}
				catch (Exception e)
				{
					System.out.println("Error parsing cigar");
					e.printStackTrace();;
					System.out.println(mate1.getStart() + " " + mate2.getStart()+" " +i+" " +mate1.getName()+" " +mate1.getaRead()); 
				}
				//full = full + " ";	
			}
		}
		else
		{
			;
		}
		if (!full.contains("~") && !full.contains("-"))
			for (int i = 0; (3*i+3) < full.length(); i++)
				countCodon(full.substring(3*i, 3*i+3), i);
	}
	
	private void writeCodons(String filename, int start, int stop)
	{
		int sum = 0; 
		try
		{
			StringTokenizer outfile = new StringTokenizer(filename, ".");
			String of = outfile.nextToken() + "_"+start+"_"+stop+"_codons.tab";
			FileWriter fw = new FileWriter(new File(of));
			
			for (int i=0; i <CODONS.length; i++)
				fw.write("\t"+CODONS[i]);
			fw.write("\tSUM\n");
			
			for(int i =0; i < codoncount.length;i++)
			{
				 fw.write(i+"\t");
				 for(int j =0; j <codoncount[i].length; j++)
				 {
					fw.write(codoncount[i][j]+"\t");
					sum = sum + codoncount[i][j];
				 }
				 fw.write(sum+"\n");
				// System.out.println(sum);
				 sum = 0;
			}
			
			fw.close();
		}
		catch (Exception e)
		{
			System.out.println("Boo");
		}
		//for (int i = 0; i < 21; i++)
			//System.out.println(AACIDS[i]+"\t"+(aacount[j][i]/sum));	
		
	}
	
	private void writeAA(String filename,int start, int stop)
	{
		int sum = 0;
		try
		{
			StringTokenizer outfile = new StringTokenizer(filename, ".");
			String of = outfile.nextToken() + "_"+start+"_"+stop+"_aa.tab";
			FileWriter fw = new FileWriter(new File(of));
			
			for (int i=0; i <AACIDS.length; i++)
				fw.write("\t"+AACIDS[i]);
			
			fw.write("\tSUM\n");
			
			for(int i =0; i < aacount.length;i++)
			{
				 fw.write(i+"\t");
				 for(int j =0; j <aacount[i].length; j++)
				 {
					fw.write(aacount[i][j]+"\t");
					sum = sum + aacount[i][j];
				 }
				fw.write(sum+"\n");
				sum = 0;
			}
			
			fw.close();
		}
		catch (Exception e)
		{
			System.out.println("Boo");
		}
		//for (int i = 0; i < 21; i++)
			//System.out.println(AACIDS[i]+"\t"+(aacount[j][i]/sum));	
		
	}
	
	private void countCodon(String codon, int i)
	{
		if (codon.equals("AAA")) //K
		{
			codoncount[i][0]++;
			aacount[i][0]++;
		}
		else if (codon.equals("AAG")) //K
		{
			codoncount[i][1]++;
			aacount[i][0]++;
		}
		else if (codon.equals("AAC")) //N
		{
			codoncount[i][2]++;
			aacount[i][1]++;
		}
		else if (codon.equals("AAT")) //N
		{
			codoncount[i][3]++;
			aacount[i][1]++;
		}
		else if (codon.equals("AGA")) //R
		{
			codoncount[i][4]++;
			aacount[i][2]++;
		}
		else if (codon.equals("AGG")) //R
		{
			codoncount[i][5]++;
			aacount[i][2]++;
		}
		else if (codon.equals("AGC")) //S
		{
			codoncount[i][6]++;
			aacount[i][3]++;
		}
		else if (codon.equals("AGT")) //S
		{
			codoncount[i][7]++;
			aacount[i][3]++;
		}
		else if (codon.equals("ACA"))//T
		{
			codoncount[i][8]++;
			aacount[i][4]++;
		}
		else if (codon.equals("ACG")) //T
		{
			codoncount[i][9]++;
			aacount[i][4]++;
		}
		else if (codon.equals("ACC")) //T
		{
			codoncount[i][10]++;
			aacount[i][4]++;
		}
		else if (codon.equals("ACT")) //T
		{
			codoncount[i][11]++;
			aacount[i][4]++;
		}
		else if (codon.equals("ATA")) //I
		{
			codoncount[i][12]++;
			aacount[i][5]++;
		}
		else if (codon.equals("ATG")) //M
		{
			codoncount[i][13]++;
			aacount[i][6]++;
		}
		else if (codon.equals("ATC")) //I
		{
			codoncount[i][14]++;
			aacount[i][5]++;
		}
		else if (codon.equals("ATT")) //I
		{
			codoncount[i][15]++;
			aacount[i][5]++;
		}
		else if (codon.equals("GAA")) //E
		{
			codoncount[i][16]++;
			aacount[i][7]++;
		}
		else if (codon.equals("GAG")) //E
		{
			codoncount[i][17]++;
			aacount[i][7]++;
		}
		else if (codon.equals("GAC"))//D
		{
			codoncount[i][18]++;
			aacount[i][8]++;
		}
		else if (codon.equals("GAT")) //D
		{
			codoncount[i][19]++;
			aacount[i][8]++;
		}
		else if (codon.equals("GGA")) //G
		{
			codoncount[i][20]++;
			aacount[i][9]++;
		}
		else if (codon.equals("GGG")) //G
		{
			codoncount[i][21]++;
			aacount[i][9]++;
		}
		else if (codon.equals("GGC")) //G
		{
			codoncount[i][22]++;
			aacount[i][9]++;
		}
		else if (codon.equals("GGT")) //G
		{
			codoncount[i][23]++;
			aacount[i][9]++;
		}
		else if (codon.equals("GCA")) //A
		{
			codoncount[i][24]++;
			aacount[i][10]++;
		}
		else if (codon.equals("GCG")) //A
		{
			codoncount[i][25]++;
			aacount[i][10]++;
		}
		else if (codon.equals("GCC")) //A
		{
			codoncount[i][26]++;
			aacount[i][10]++;
		}
		else if (codon.equals("GCT")) //A
		{
			codoncount[i][27]++;
			aacount[i][10]++;
		}
		else if (codon.equals("GTA")) //V
		{
			codoncount[i][28]++;
			aacount[i][11]++;
		}
		else if (codon.equals("GTG")) //V
		{
			codoncount[i][29]++;
			aacount[i][11]++;
		}
		else if (codon.equals("GTC")) //V
		{
			codoncount[i][30]++;
			aacount[i][11]++;
		}
		else if (codon.equals("GTT")) //V
		{
			codoncount[i][31]++;
			aacount[i][11]++;
		}
		else if (codon.equals("CAA")) //Q
		{
			codoncount[i][32]++;
			aacount[i][12]++;
		}
		else if (codon.equals("CAG")) //Q
		{
			codoncount[i][33]++;
			aacount[i][12]++;
		}
		else if (codon.equals("CAC")) //H
		{
			codoncount[i][34]++;
			aacount[i][13]++;
		}
		else if (codon.equals("CAT")) //H
		{
			codoncount[i][35]++;
			aacount[i][13]++;
		}
		else if (codon.equals("CGA")) //R
		{
			codoncount[i][36]++;
			aacount[i][2]++;
		}
		else if (codon.equals("CGG")) //R
		{
			codoncount[i][37]++;
			aacount[i][2]++;
		}
		else if (codon.equals("CGC")) //R
		{
			codoncount[i][38]++;
			aacount[i][2]++;
		}
		else if (codon.equals("CGT")) //R
		{
			codoncount[i][39]++;
			aacount[i][2]++;
		}else if (codon.equals("CCA")) //P
		{
			codoncount[i][40]++;
			aacount[i][14]++;
		}
		else if (codon.equals("CCG")) //P
		{
			codoncount[i][41]++;
			aacount[i][14]++;
		}
		else if (codon.equals("CCC")) //P
		{
			codoncount[i][42]++;
			aacount[i][14]++;
		}
		else if (codon.equals("CCT")) //P
		{
			codoncount[i][43]++;
			aacount[i][14]++;
		}
		else if (codon.equals("CTA")) //L
		{
			codoncount[i][44]++;
			aacount[i][15]++;
		}
		else if (codon.equals("CTG")) //L
		{
			codoncount[i][45]++;
			aacount[i][15]++;
		}
		else if (codon.equals("CTC")) //L
		{
			codoncount[i][46]++;
			aacount[i][15]++;
		}
		else if (codon.equals("CTT")) //L
		{
			codoncount[i][47]++;
			aacount[i][15]++;
		}
		else if (codon.equals("TAA")) //*
		{
			codoncount[i][48]++;
			aacount[i][16]++;
		}
		else if (codon.equals("TAG")) //*
		{
			codoncount[i][49]++;
			aacount[i][16]++;
		}
		else if (codon.equals("TAC")) //Y
		{
			codoncount[i][50]++;
			aacount[i][17]++;
		}
		else if (codon.equals("TAT")) //Y
		{
			codoncount[i][51]++;
			aacount[i][17]++;
		}
		else if (codon.equals("TGA")) //*
		{
			codoncount[i][52]++;
			aacount[i][16]++;
		}
		else if (codon.equals("TGG")) //W
		{
			codoncount[i][53]++;
			aacount[i][18]++;
		}
		else if (codon.equals("TGC")) //C
		{
			codoncount[i][54]++;
			aacount[i][19]++;
		}
		else if (codon.equals("TGT")) //C
		{
			codoncount[i][55]++;
			aacount[i][19]++;
		}
		else if (codon.equals("TCA")) //S
		{
			codoncount[i][56]++;
			aacount[i][3]++;
		}
		else if (codon.equals("TCG")) //S
		{
			codoncount[i][57]++;
			aacount[i][3]++;
		}
		else if (codon.equals("TCC")) //S
		{
			codoncount[i][58]++;
			aacount[i][3]++;
		}
		else if (codon.equals("TCT")) //S
		{
			codoncount[i][59]++;
			aacount[i][3]++;
		}
		else if (codon.equals("TTA")) //L
		{
			codoncount[i][60]++;
			aacount[i][15]++;
		}
		else if (codon.equals("TTG")) //L
		{
			codoncount[i][61]++;
			aacount[i][15]++;
		}
		else if (codon.equals("TTC")) //F
		{
			codoncount[i][62]++;
			aacount[i][20]++;
		}
		else if (codon.equals("TTT")) //F
		{
			codoncount[i][63]++;
			aacount[i][20]++;
		}	
	}

	private class SAMString
	{
		private String name; //read name (cluster ID)
		private String flag; //flag
		private String gname;
		private int start; //mapping start
		private int mq; //mapping quality
		private String cigar; //cigar mutations from ref
		private String mate; //name of mate
		private int mpos; //mpos mapping
		private int tlen; //frag length
		private String read; //actual read
		private String qual; //actual quality 
		private String aligned_read;
		
		public SAMString(String n, String f, String g, int s, int mpq, String cig, String mname, int mp, int t, String r, String q)
		{
			name = n;
			flag = f;
			gname = g;
			start = s;
			mq = mpq;
			cigar = cig;
			mate = mname;
			mpos = mp;
			tlen = t;
			read = r;
			qual = q;
			aligned_read = parseCIGAR();
		}
		
		public String getaRead()
		{
			return aligned_read;
		}
		
		public String parseCIGAR()
		{
			String m1 = read;
			
			StringTokenizer ct = new StringTokenizer(cigar, "DHIMNPSX=", true);
			
			int j = 0;
			
			while (ct.hasMoreTokens())
			{
				String tmp = ct.nextToken();
				
				if (!tmp.equals("*"))
				{
					int i = Integer.valueOf(tmp); //cigar int
					int k=j+i-1; //position in string, cigar is +1 instead of zero
					
					String delim = ct.nextToken();
		
					if (delim.equals("D"))
					{
						String s = "";
						for (int a=0; a<i; a++)
						{
							s = s + "-";
						}
						m1 = m1.substring(0,k) + s + m1.substring(k+1,m1.length());
					}
					else if (delim.equals("S"))
					{
						String s = "";
						for (int a=0; a<i; a++)
						{
							s = s + "N";
						}
						m1 = m1.substring(0,j) + s + m1.substring(k+1,m1.length());
					}
					else if (delim.equals("I"))
					{
						m1 = m1.substring(0,j-1) + "~" + m1.substring(k, m1.length());
					}
					
					j = k;
				}
			}			
			m1 = m1.replace("N","");
			return m1;
		}
		
		
		public String getName()
		{
			return name;	
		}
		
		public String getRead()
		{
			return read;
		}
		
		public String getQual()
		{
			return qual;
		}
		
		public int getTLength()
		{
			return tlen;
		}
		
		public int getStart()
		{
			return start;
		}
		
		public int getEnd()
		{
			return mpos;
		}
		
		private String reverseComplementRead()
		{
			String b = new StringBuilder(read).reverse().toString();
			b = b.replace("A", "W");
			b = b.replace("T", "X");
			b = b.replace("G", "Y");
			b = b.replace("C", "Z");
			b = b.replace("W", "T");
			b = b.replace("X", "A");
			b = b.replace("Y", "C");
			b = b.replace("Z", "G");
			
			return b;
		}
		
		public String getCig()
		{
			return cigar;
		}
	}
	
	private class SAMJoin
	{
		private SAMString mate1;
		private SAMString mate2;
		private int mstart;
		private int mend;
		private String name;
		
		public SAMJoin(SAMString m1, SAMString m2)
		{
			mate1 = m1;
			mate2 = m2;
			mstart = m1.getStart();
			mend = m1.getEnd();
			name = m1.getName();
		}
		
		public SAMString getMate1()
		{
			return mate1;
		}
		
		public SAMString getMate2()
		{
			return mate2;
		}
		
		public String getName()
		{
			return name;
		}
		
		public int getEnd()
		{
			return (mstart + mate1.getTLength());
		}
		
		public int getStart()
		{
			return mstart;
		}
	}
}