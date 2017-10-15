/*
 * Program name: Sequence.java
 * Date: 2016/09/30
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		Sequence for [dNTP1/dNTP2], [-/dNTP1/.../dNTPn] or IUPAC format.
 */

package bio.rehunt.seq;

import java.util.*;

/**
 * Sequence for [dNTP1/dNTP2], [-/dNTP1/.../dNTPn] or IUPAC format.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 */
public class Sequence {
	private String iupac_seq = null;
	private String dntp_seq = null;
	private LinkedList<String> iupac_list = null;
	private LinkedList<String> dntps_list = null;
	private LinkedList<Integer> var_pos_list = null;

	/**
	 * Judge if sequence contains dNTPs or IUPAC format.
	 * @param seq Seqeuence.
	 * @return is_var_seq True is variation sequence, else false.
	 */
	public boolean isVarSeq(String seq) {
		boolean is_var_seq = false;
		String upper_seq = seq.toUpperCase();
		StringBuffer strBuff_seq = new StringBuffer();
		LinkedList<String[]> list_sym = new LinkedList<String[]>();
		for(int i=0;i<upper_seq.length();i++) {
			if(upper_seq.charAt(i)=='A' || upper_seq.charAt(i)=='T' || upper_seq.charAt(i)=='C' || upper_seq.charAt(i)=='G')
				strBuff_seq.append(upper_seq.charAt(i));
			else if(upper_seq.charAt(i)=='M' || upper_seq.charAt(i)=='R' || upper_seq.charAt(i)=='W' || upper_seq.charAt(i)=='S' ||
				upper_seq.charAt(i)=='Y' || upper_seq.charAt(i)=='K' || upper_seq.charAt(i)=='V' || upper_seq.charAt(i)=='H' ||
				upper_seq.charAt(i)=='D' || upper_seq.charAt(i)=='B' || upper_seq.charAt(i)=='N' ||
				upper_seq.charAt(i)=='[' || upper_seq.charAt(i)=='/' || upper_seq.charAt(i)==']' || upper_seq.charAt(i)=='-') {
				if(upper_seq.charAt(i)=='[' || upper_seq.charAt(i)==']') {
					String[] sym = new String[2];
					sym[0] = Character.toString(upper_seq.charAt(i));
					sym[1] = Integer.toString(strBuff_seq.length());
					list_sym.add(sym);
				}
				strBuff_seq.append(upper_seq.charAt(i));
				is_var_seq = true;
			}
		}
		// check the symbol '[' and ']'
		//int leftSym = -1;	// the position of '['
		//int rightSym = -1;	// the position of ']'
		boolean check_ok = false;
		if(list_sym.size() != 0) {
			for(int i=0;i<list_sym.size();i++) {
				String[] sym = list_sym.get(i);
				// symbol '['
				if(i%2 == 0 && sym[0].equals("[")) {
					//leftSym = Integer.parseInt(sym[1]);
					//System.out.println(leftSym);
					continue;
				}
				// symbol ']'
				else if(i%2 == 1 && sym[0].equals("]")) {
					//rightSym = Integer.parseInt(sym[1]);
					//System.out.println(rightSym);
					check_ok = true;
				}
				else {
					System.out.println("Symbol '[' or ']' error!");
					is_var_seq = false;
					return is_var_seq;
				}
				if(check_ok) {	// check symbol '[' and ']' ok
					if(i%2 == 1 && i != list_sym.size()-1)
						check_ok = false;
				}
			}
			if(!check_ok) {
				System.out.println("Symbol '[' or ']' error!");
				is_var_seq = false;
			}
		}
		return is_var_seq;
	}

	/**
	 * Make sequence with [dNTP1/dNTP2], [-/dNTP1/.../dNTPn] or IUPAC format.
	 * @param seq The seqeuence with IUPAC format or [dNTP1/dNTP1/.../dNTPn] format.
	 * @return If make sequence successfully return true, else return false.
	 */
	public boolean makeVarSeq(String seq) {
		String upper_seq = seq.toUpperCase();
		// trans IUPAC to [dNTP1/dNTP2]
		StringBuffer strBuff_dNTPseq = new StringBuffer();
		LinkedList<String[]> list_sym = new LinkedList<String[]>();
		for(int i=0;i<upper_seq.length();i++) {
			if(upper_seq.charAt(i)=='A' || upper_seq.charAt(i)=='T' || upper_seq.charAt(i)=='C' || upper_seq.charAt(i)=='G')
				strBuff_dNTPseq.append(upper_seq.charAt(i));
			if(upper_seq.charAt(i)=='M' || upper_seq.charAt(i)=='R' || upper_seq.charAt(i)=='W' || upper_seq.charAt(i)=='S' ||
				upper_seq.charAt(i)=='Y' || upper_seq.charAt(i)=='K' || upper_seq.charAt(i)=='V' || upper_seq.charAt(i)=='H' ||
				upper_seq.charAt(i)=='D' || upper_seq.charAt(i)=='B' || upper_seq.charAt(i)=='N') {
				String dNTPs = IUPACTodNTPs(Character.toString(upper_seq.charAt(i)));
				strBuff_dNTPseq.append(dNTPs);
			}
			else if(upper_seq.charAt(i)=='[' || upper_seq.charAt(i)=='/' || upper_seq.charAt(i)==']' || upper_seq.charAt(i)=='-')
				strBuff_dNTPseq.append(upper_seq.charAt(i));
		}
		dntp_seq = strBuff_dNTPseq.toString();
		// filter sequence
		StringBuffer strBuff_seq = new StringBuffer();
		for(int i=0;i<dntp_seq.length();i++) {
			if( dntp_seq.charAt(i)=='A' || dntp_seq.charAt(i)=='T' || dntp_seq.charAt(i)=='C' || dntp_seq.charAt(i)=='G' ||
				dntp_seq.charAt(i)=='M' || dntp_seq.charAt(i)=='R' || dntp_seq.charAt(i)=='W' || dntp_seq.charAt(i)=='S' ||
				dntp_seq.charAt(i)=='Y' || dntp_seq.charAt(i)=='K' || dntp_seq.charAt(i)=='V' || dntp_seq.charAt(i)=='H' ||
				dntp_seq.charAt(i)=='D' || dntp_seq.charAt(i)=='B' || dntp_seq.charAt(i)=='N' ||
				dntp_seq.charAt(i)=='[' || dntp_seq.charAt(i)=='/' || dntp_seq.charAt(i)==']' || dntp_seq.charAt(i)=='-') {
				if(dntp_seq.charAt(i)=='[' || dntp_seq.charAt(i)==']') {
					String[] sym = new String[2];
					sym[0] = Character.toString(dntp_seq.charAt(i));
					sym[1] = Integer.toString(strBuff_seq.length());
					list_sym.add(sym);
				}
				strBuff_seq.append(dntp_seq.charAt(i));
			}
		}
		iupac_seq = strBuff_seq.toString();
		// check the symbol '[' and ']'
		iupac_list = new LinkedList<String>();
		dntps_list = new LinkedList<String>();
		var_pos_list = new LinkedList<Integer>();
		StringBuffer strBuff_tmp = new StringBuffer(iupac_seq);
		int leftSym = -1;	// the position of '['
		int rightSym = -1;	// the position of ']'
		boolean check_ok = false;
		if(list_sym.size() != 0) {
			for(int i=0;i<list_sym.size();i++) {
				String[] sym = list_sym.get(i);
				// symbol '['
				if(i%2 == 0 && sym[0].equals("[")) {
					leftSym = Integer.parseInt(sym[1]);
					//System.out.println(leftSym);
				}
				// symbol ']'
				else if(i%2 == 1 && sym[0].equals("]")) {
					rightSym = Integer.parseInt(sym[1]);
					//System.out.println(rightSym);
					check_ok = true;
				}
				else {
					System.out.println("Symbol '[' or ']' error!");
					return false;
				}
				if(check_ok) {	// check symbol '[' and ']' ok
					// replace [dNTP1/dNTP2] to IUPAC
					StringBuffer strBuff_iupac = new StringBuffer();
					String dNTPs = iupac_seq.substring(leftSym, rightSym+1);
					String iupac = dNTPsToIUPAC(dNTPs);
					dntps_list.add(dNTPs);
					iupac_list.add(iupac);
					for(int j=0;j<dNTPs.length();j++)
						strBuff_iupac.append(" ");
					strBuff_iupac.replace(0, 1, iupac);
					strBuff_tmp.replace(leftSym, rightSym+1, strBuff_iupac.toString());
					if(i%2 == 1 && i != list_sym.size()-1)
						check_ok = false;
				}
			}
			if(!check_ok) {
				System.out.println("Symbol '[' or ']' error!");
				return false;
			}
			String strTmp = strBuff_tmp.toString();
			StringBuffer strBuff_ok = new StringBuffer();
			int spaceNum = 0;
			for(int i=0;i<strTmp.length();i++) {
				if(strTmp.charAt(i) != ' ') {
					strBuff_ok.append(strTmp.charAt(i));
					if(strTmp.charAt(i)=='M' || strTmp.charAt(i)=='R' || strTmp.charAt(i)=='W' || strTmp.charAt(i)=='S' ||
						strTmp.charAt(i)=='Y' || strTmp.charAt(i)=='K' || strTmp.charAt(i)=='V' || strTmp.charAt(i)=='H' ||
						strTmp.charAt(i)=='D' || strTmp.charAt(i)=='B' || strTmp.charAt(i)=='N') {
						int snp_pos = i - spaceNum;
						var_pos_list.add(snp_pos);
					}
				}
				else
					spaceNum++;
			}
			iupac_seq = strBuff_ok.toString();
		}
		return true;
	}

	/**
	 * Get allele bases from dNTPs format.
	 * @param dNTPs dNTPs format.
	 * @return Allele bases array.
	 */
	public String[] getAlleleSplit(String dNTPs) {
		String[] allele_split = dNTPs.split("/");
		allele_split[0] = allele_split[0].substring(1);
		allele_split[allele_split.length-1] = allele_split[allele_split.length-1].substring(0, allele_split[allele_split.length-1].length() - 1);
		return allele_split;
	}

	/**
	 * Transform sequence from dNTPs format to IUPAC foramt.
	 * @param dNTPs dNTPs format.
	 * @return IUPAC foramt.
	 */
	public String dNTPsToIUPAC(String dNTPs) {
		String IUPAC = "";
		if(dNTPs.equals("[A/C]") || dNTPs.equals("[C/A]"))
			IUPAC = "M";
		else if(dNTPs.equals("[A/G]") || dNTPs.equals("[G/A]"))
			IUPAC = "R";
		else if(dNTPs.equals("[A/T]") || dNTPs.equals("[T/A]"))
			IUPAC = "W";
		else if(dNTPs.equals("[C/G]") || dNTPs.equals("[G/C]"))
			IUPAC = "S";
		else if(dNTPs.equals("[C/T]") || dNTPs.equals("[T/C]"))
			IUPAC = "Y";
		else if(dNTPs.equals("[G/T]") || dNTPs.equals("[T/G]"))
			IUPAC = "K";
		else if(dNTPs.equals("[A/C/G]") || dNTPs.equals("[A/C/C]") ||
		        dNTPs.equals("[C/A/G]") || dNTPs.equals("[C/G/A]") ||
		        dNTPs.equals("[G/A/C]") || dNTPs.equals("[G/C/A]"))
			IUPAC = "V";
		else if(dNTPs.equals("[A/C/T]") || dNTPs.equals("[A/T/C]") ||
		        dNTPs.equals("[C/A/T]") || dNTPs.equals("[C/T/A]") ||
		        dNTPs.equals("[G/A/T]") || dNTPs.equals("[G/T/A]"))
			IUPAC = "H";
		else if(dNTPs.equals("[A/G/T]") || dNTPs.equals("[A/T/G]") ||
		        dNTPs.equals("[G/A/T]") || dNTPs.equals("[G/T/A]") ||
		        dNTPs.equals("[T/A/G]") || dNTPs.equals("[T/G/A]"))
			IUPAC = "D";
		else if(dNTPs.equals("[C/G/T]") || dNTPs.equals("[C/T/G]") ||
		        dNTPs.equals("[G/C/T]") || dNTPs.equals("[G/T/C]") ||
		        dNTPs.equals("[T/C/G]") || dNTPs.equals("[T/G/C]"))
			IUPAC = "B";
		else
			IUPAC = "N";

		return IUPAC;
	}

	/**
	 * Transform sequence from IUPAC format to dNTPs foramt.
	 * @param IUPAC IUPAC format.
	 * @return dNTPs format.
	 */
	public String IUPACTodNTPs(String IUPAC) {
		String dNTPs = "";
		if(IUPAC.equals("M"))
			dNTPs = "[A/C]";
		else if(IUPAC.equals("R"))
			dNTPs = "[A/G]";
		else if(IUPAC.equals("W"))
			dNTPs = "[A/T]";
		else if(IUPAC.equals("S"))
			dNTPs = "[C/G]";
		else if(IUPAC.equals("Y"))
			dNTPs = "[C/T]";
		else if(IUPAC.equals("K"))
			dNTPs = "[G/T]";
		else if(IUPAC.equals("V"))
			dNTPs = "[A/C/G]";
		else if(IUPAC.equals("H"))
			dNTPs = "[A/C/T]";
		else if(IUPAC.equals("D"))
			dNTPs = "[A/G/T]";
		else if(IUPAC.equals("B"))
			dNTPs = "[C/G/T]";
		else if(IUPAC.equals("N"))
			dNTPs = "[A/C/G/T]";
		return dNTPs;
	}

	/**
	 * Transform sequence into a reverse sequence.
	 * @param seq The seqeuence.
	 * @return The reverse sequence.
	 */
	public String reverseTrans(String seq) {
		seq = seq.toUpperCase();
		StringBuffer strBuff_seq = new StringBuffer(seq);
		return strBuff_seq.reverse().toString();
	}

	/**
	 * Transform sequence into a complementary sequence.
	 * @param seq IUPAC format sequence.
	 * @return The complementary sequence.
	 */
	public String complementaryTrans(String seq) {
		seq = seq.toUpperCase();
		StringBuffer strBuff_complementarySeq = new StringBuffer();
		for(int i=0;i<seq.length();i++) {
			if(seq.charAt(i) == 'A')
				strBuff_complementarySeq.append("T");
			else if(seq.charAt(i) == 'C')
				strBuff_complementarySeq.append("G");
			else if(seq.charAt(i) == 'G')
				strBuff_complementarySeq.append("C");
			else if(seq.charAt(i) == 'T')
				strBuff_complementarySeq.append("A");
			else if(seq.charAt(i) == 'M')
				strBuff_complementarySeq.append("K");
			else if(seq.charAt(i) == 'R')
				strBuff_complementarySeq.append("Y");
			else if(seq.charAt(i) == 'W')
				strBuff_complementarySeq.append("W");
			else if(seq.charAt(i) == 'S')
				strBuff_complementarySeq.append("S");
			else if(seq.charAt(i) == 'Y')
				strBuff_complementarySeq.append("R");
			else if(seq.charAt(i) == 'K')
				strBuff_complementarySeq.append("M");
			else if(seq.charAt(i) == 'V')
				strBuff_complementarySeq.append("B");
			else if(seq.charAt(i) == 'H')
				strBuff_complementarySeq.append("D");
			else if(seq.charAt(i) == 'D')
				strBuff_complementarySeq.append("H");
			else if(seq.charAt(i) == 'B')
				strBuff_complementarySeq.append("V");
			else if(seq.charAt(i) == 'N')
				strBuff_complementarySeq.append("N");
			else if(seq.charAt(i) == '[')
				strBuff_complementarySeq.append("[");
			else if(seq.charAt(i) == '/')
				strBuff_complementarySeq.append("/");
			else if(seq.charAt(i) == ']')
				strBuff_complementarySeq.append("]");
			else if(seq.charAt(i) == '-')
				strBuff_complementarySeq.append("-");
		}
		return strBuff_complementarySeq.toString();
	}

	/**
	 * Transform sequence into an anti-sense sequence (reverse and complementary sequence).
	 * @param seq IUPAC format sequence.
	 * @return The reverse and complementary sequence.
	 */
	public String antiTrans(String seq) {
		seq = seq.toUpperCase();
		StringBuffer strBuff_antiseq = new StringBuffer();
		for(int i=seq.length()-1;i>=0;i--) {
			if(seq.charAt(i) == 'A')
				strBuff_antiseq.append("T");
			else if(seq.charAt(i) == 'C')
				strBuff_antiseq.append("G");
			else if(seq.charAt(i) == 'G')
				strBuff_antiseq.append("C");
			else if(seq.charAt(i) == 'T')
				strBuff_antiseq.append("A");
			else if(seq.charAt(i) == 'M')
				strBuff_antiseq.append("K");
			else if(seq.charAt(i) == 'R')
				strBuff_antiseq.append("Y");
			else if(seq.charAt(i) == 'W')
				strBuff_antiseq.append("W");
			else if(seq.charAt(i) == 'S')
				strBuff_antiseq.append("S");
			else if(seq.charAt(i) == 'Y')
				strBuff_antiseq.append("R");
			else if(seq.charAt(i) == 'K')
				strBuff_antiseq.append("M");
			else if(seq.charAt(i) == 'V')
				strBuff_antiseq.append("B");
			else if(seq.charAt(i) == 'H')
				strBuff_antiseq.append("D");
			else if(seq.charAt(i) == 'D')
				strBuff_antiseq.append("H");
			else if(seq.charAt(i) == 'B')
				strBuff_antiseq.append("V");
			else if(seq.charAt(i) == 'N')
				strBuff_antiseq.append("N");
			else if(seq.charAt(i) == ']')
				strBuff_antiseq.append("[");
			else if(seq.charAt(i) == '/')
				strBuff_antiseq.append("/");
			else if(seq.charAt(i) == '[')
				strBuff_antiseq.append("]");
			else if(seq.charAt(i) == '-')
				strBuff_antiseq.append("-");
		}
		return strBuff_antiseq.toString();
	}

	/**
	 * Get IUPAC format sequence.
	 * @return IUPAC format sequence.
	 */
	public String getIUPACSeq() {
		return iupac_seq;
	}

	/**
	 * Get dNTPs format sequence.
	 * @return dNTPs format sequence.
	 */
	public String getdNTPSeq() {
		return dntp_seq;
	}

	/**
	 * Get IUPAC list after the method makeVarSeq(String seq).
	 * @return IUPAC list in the sequence.
	 */
	// get iupac list
	public List<String> getIUPACList() {
		return iupac_list;
	}

	/**
	 * Get dNTPs list after the method makeVarSeq(String seq).
	 * @return dNTPs list in the sequence.
	 */
	public List<String> getdNTPsList() {
		return dntps_list;
	}

	/**
	 * Get variation position list after the method makeVarSeq(String seq).
	 * @return variation position list in the sequence.
	 */
	public List<Integer> getVarPosList() {
		return var_pos_list;
	}

/*	public static void main(String args[]) throws Exception {
		Sequence sequence = new Sequence();
		String seq = "ANARC[-/G]TCAT[-/GC]TGCA";
		if(sequence.isVarSeq(seq)) {
			boolean makeSuccess = sequence.makeVarSeq(seq);
			if(makeSuccess) {
				String iupac_seq = sequence.getIUPACSeq();
				String dntp_seq = sequence.getdNTPSeq();
				List<String> iupac_list = sequence.getIUPACList();
				List<String> dntps_list = sequence.getdNTPsList();
				List<Integer> var_pos_list = sequence.getVarPosList();
				System.out.println("IUPAC SNP: " + iupac_seq);
				System.out.println("dNTPs SNP: " + dntp_seq);
				for(int i=0;i<dntps_list.size();i++) {
					System.out.println(iupac_list.get(i));
					System.out.println(dntps_list.get(i));
					System.out.println(var_pos_list.get(i));
					String dntps = dntps_list.get(i);
					String[] allele_split = sequence.getAlleleSplit(dntps);
					for(int j=0;j<allele_split.length;j++)
						System.out.println("allele_split: " + allele_split[j]);
				}
			}
		}
	}*/
}