/*
 * Program name: SeqProcess.java
 * Date: 2016/09/30
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		Sequence process function.
 */

package bio.rehunt.seq;

import java.util.*;

/**
 * Sequence process function.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 */
public class SeqProcess {
	/**
	 * Get the first IUPAC position from a sequence.
	 * @param seq Sequence with IUPAC format.
	 * @return The first IUPAC position in the sequence. If no IUPAC in the sequence then return -1.
	 */
	public int getFirstIUPACPos(String seq) {
		int iupac_pos = -1;	// the IUPAC position
		for(int i=0;i<seq.length();i++) {
			if(seq.charAt(i)=='M' || seq.charAt(i)=='R' || seq.charAt(i)=='W' || seq.charAt(i)=='S' ||
				seq.charAt(i)=='Y' || seq.charAt(i)=='K' || seq.charAt(i)=='V' || seq.charAt(i)=='H' ||
				seq.charAt(i)=='D' || seq.charAt(i)=='B' || seq.charAt(i)=='N') {
				iupac_pos = i;
				break;
			}
		}
		return iupac_pos;
	}

	/**
	 * Get flanking sequence.
	 * @param seq Sequence with IUPAC format.
	 * @param flank_len Flanking sequence length.
	 * @param central_pos The central position.
	 * @return An string array for a flanking sequence (array[0]) and its central position (array[1]).
	 */
	public String[] getFlankSeq(String seq, int flank_len, int central_pos) {
		int central_pos_flank = -1;
		StringBuffer strBuff_seqFlank = new StringBuffer();
		int len = flank_len;
		int pos = central_pos;
		if(pos == -1)	// if pos=-1, set central_pos be the first IUPAC position
			pos = getFirstIUPACPos(seq);
		String seq5 = "";
		String seq3 = "";
		if(pos-len >= 0) {
			seq5 = seq.substring(pos-len, pos);
			central_pos_flank = len;
		}
		else {
			seq5 = seq.substring(0, pos);
			central_pos_flank = pos;
		}
		if(pos+len < seq.length())
			seq3 = seq.substring(pos+1, pos+len+1);
		else
			seq3 = seq.substring(pos+1);
		strBuff_seqFlank.append(seq5);
		strBuff_seqFlank.append(seq.charAt(pos));
		strBuff_seqFlank.append(seq3);
		String[] return_data = new String[2];
		return_data[0] = strBuff_seqFlank.toString();	// flank sequence
		return_data[1] = Integer.toString(central_pos_flank);	// central position in the flank sequence
		return return_data;
	}

	/**
	 * Get the multiple sequences according to all dNTPs.
	 * @param seq Sequence with IUPAC format.
	 * @param dntps_list dNTPs list in the sequence.
	 * @param dntps_index The index of dNTPs list, the dntps_index must be set as 0.
	 * @return A sequence list according to all dNTPs.
	 */
	public List<String> getMultiSeq(String seq, List<String> dntps_list, int dntps_index) {
		List<String> seq_list = new LinkedList<String>();
		List<String> seqs = null;
		for(int i=0;i<seq.length();i++) {
			if(seq.charAt(i)=='M' || seq.charAt(i)=='R' || seq.charAt(i)=='W' || seq.charAt(i)=='S' ||
				seq.charAt(i)=='Y' || seq.charAt(i)=='K' || seq.charAt(i)=='V' || seq.charAt(i)=='H' ||
				seq.charAt(i)=='D' || seq.charAt(i)=='B' || seq.charAt(i)=='N') {
				if(dntps_index < dntps_list.size()) {
					String dntps = dntps_list.get(dntps_index);
					String allele = "";
					int current_count = 0;
					for(int j=0;j<dntps.length();j++) {
						if(dntps.charAt(j) == '[')
							current_count = j;
						else if(dntps.charAt(j) == '/' || dntps.charAt(j) == ']') {
							allele = dntps.substring(current_count+1, j);
							StringBuffer strBuff = new StringBuffer(seq);
							seqs = getMultiSeq(strBuff.replace(i, i+1, allele).toString(), dntps_list, dntps_index+1);
							for(int k=0;k<seqs.size();k++)
								seq_list.add(seqs.get(k).toString());
							current_count = j;
						}
					}
				}
				break;
			}
		}
		if(seq_list.size() == 0)
			seq_list.add(seq);
		return seq_list;
	}

	/**
	 * Get the multiple sequences according to the assigned position.
	 * @param seq Sequence with IUPAC format.
	 * @param pos The assigned position.
	 * @return Sequence list according to the assigned position.
	 */
	public List<String> getPosMultiSeq(String seq, int pos) {
		List<String> seq_list = new LinkedList<String>();

		StringBuffer iupacSNP_strBuff = new StringBuffer(seq);
		if(seq.charAt(pos)=='M') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "A").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "C").toString();
			seq_list.add(seq2);
		}
		else if(seq.charAt(pos)=='R') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "A").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "G").toString();
			seq_list.add(seq2);
		}
		else if(seq.charAt(pos)=='W') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "A").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "T").toString();
			seq_list.add(seq2);
		}
		else if(seq.charAt(pos)=='S') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "C").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "G").toString();
			seq_list.add(seq2);
		}
		else if(seq.charAt(pos)=='Y') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "C").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "T").toString();
			seq_list.add(seq2);
		}
		else if(seq.charAt(pos)=='K') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "G").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "T").toString();
			seq_list.add(seq2);
		}
		else if(seq.charAt(pos)=='V') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "A").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "C").toString();
			seq_list.add(seq2);
			String seq3 = iupacSNP_strBuff.replace(pos, pos+1, "G").toString();
			seq_list.add(seq3);
		}
		else if(seq.charAt(pos)=='H') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "A").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "C").toString();
			seq_list.add(seq2);
			String seq3 = iupacSNP_strBuff.replace(pos, pos+1, "T").toString();
			seq_list.add(seq3);
		}
		else if(seq.charAt(pos)=='D') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "A").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "G").toString();
			seq_list.add(seq2);
			String seq3 = iupacSNP_strBuff.replace(pos, pos+1, "T").toString();
			seq_list.add(seq3);
		}
		else if(seq.charAt(pos)=='B') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "C").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "G").toString();
			seq_list.add(seq2);
			String seq3 = iupacSNP_strBuff.replace(pos, pos+1, "T").toString();
			seq_list.add(seq3);
		}
		else if(seq.charAt(pos)=='N') {
			String seq1 = iupacSNP_strBuff.replace(pos, pos+1, "A").toString();
			seq_list.add(seq1);
			String seq2 = iupacSNP_strBuff.replace(pos, pos+1, "C").toString();
			seq_list.add(seq2);
			String seq3 = iupacSNP_strBuff.replace(pos, pos+1, "G").toString();
			seq_list.add(seq3);
			String seq4 = iupacSNP_strBuff.replace(pos, pos+1, "T").toString();
			seq_list.add(seq4);
		}
		if(seq_list.size() == 0)	// return original sequence
			seq_list.add(seq);
		return seq_list;
	}

	/**
	 * Get the multiple sequences with dNTPs format according to the assigned position and dNTPs.
	 * @param seq Sequence with IUPAC format.
	 * @param pos The assigned position in the sequence.
	 * @param dntps_list The dNTPs list in the sequence.
	 * @param var_pos_list The variation position list in the sequence.
	 * @return Sequence list with dNTPs format according to the assigned position and dNTPs.
	 */
	public List<String> getPosMultiSeq_dNTPs(String seq, int pos, List<String> dntps_list, List<Integer> var_pos_list) {
		List<String> seq_list = new LinkedList<String>();

		if(pos == -1)	// if pos=-1, set pos be the first IUPAC position
			pos = getFirstIUPACPos(seq);

		// find dNTPs from variation position list by assigned position
		String dntps = "";
		for(int i=0;i<var_pos_list.size();i++) {
			int pos_tmp = var_pos_list.get(i);
			if(pos_tmp == pos) {
				dntps = dntps_list.get(i);
				break;
			}
		}
		StringBuffer seq_strBuff = new StringBuffer(seq);
		int start_pos = 1;	// start position after symbol '['
		for(int i=0;i<dntps.length();i++) {
			if(dntps.charAt(i) == '/' || dntps.charAt(i) == ']') {
				String seq_tmp = seq_strBuff.replace(pos, pos+1, dntps.substring(start_pos, i)).toString();
				seq_list.add(seq_tmp);
				start_pos = i + 1;
			}
		}
		if(seq_list.size() == 0)	// return original sequence
			seq_list.add(seq);
		return seq_list;
	}

	/**
	 * Get the multiple sequences according to the assigned position and dNTPs.
	 * @param seq Sequence with IUPAC format.
	 * @param pos The assigned position in the sequence.
	 * @param dntps The dNTPs in the sequence.
	 * @return Sequence list according to the assigned position and dNTPs.
	 */
	public List<String> getPosMultiSeq_dNTPs(String seq, int pos, String dntps) {
		List<String> seq_list = new LinkedList<String>();
		if(pos == -1)	// if pos=-1, set pos be the first IUPAC position
			pos = getFirstIUPACPos(seq);
		StringBuffer iupacSNP_strBuff = new StringBuffer(seq);
		int start_pos = 1;	// start position after symbol '['
		for(int i=0;i<dntps.length();i++) {
			if(dntps.charAt(i) == '/' || dntps.charAt(i) == ']') {
				String seq_tmp = iupacSNP_strBuff.replace(pos, pos+1, dntps.substring(start_pos, i)).toString();
				seq_list.add(seq_tmp);
				start_pos = i + 1;
			}
		}
		if(seq_list.size() == 0)	// return original sequence
			seq_list.add(seq);
		return seq_list;
	}

/*	public static void main(String args[]) throws Exception {
		// sequence process
		Sequence sequence = new Sequence();
		String seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAWAA[-/AGC]AA";
		String iupac_seq = null;
		String dntp_seq = null;
		List<String> iupac_list = null;
		List<String> dntps_list = null;
		List<Integer> var_pos_list = null;
		if(sequence.isVarSeq(seq)) {
			boolean is_success = sequence.makeVarSeq(seq);
			if(is_success) {
				iupac_seq = sequence.getIUPACSeq();
				dntp_seq = sequence.getdNTPSeq();
				iupac_list = sequence.getIUPACList();
				dntps_list = sequence.getdNTPsList();
				var_pos_list = sequence.getVarPosList();
			}
		}
		System.out.println("iupac_seq: " + iupac_seq);
		// get multiSeq for all dNTPs
		SeqProcess seqProcess = new SeqProcess();
		int dntps_index = 0;
		List<String> multiSeq = seqProcess.getMultiSeq(iupac_seq, dntps_list, dntps_index);
		for(int i=0;i<multiSeq.size();i++)
			System.out.println("multiSeq: " + multiSeq.get(i));
		// get sequence data
		int flank_len = 10;
		int central_pos = -1;
		String[] seq_data = seqProcess.getFlankSeq(iupac_seq, flank_len, central_pos);
		String seq_flank = seq_data[0];
		int central_pos_flank = Integer.parseInt(seq_data[1]);
		System.out.println("seq_flank: " + seq_flank);
		System.out.println("central_pos_flank: " + central_pos_flank);
		// get posMultiSeq from IUPAC seq for part flanking sequence
		List<String> posMultiSeq = seqProcess.getPosMultiSeq(seq_flank, central_pos_flank);
		for(int i=0;i<posMultiSeq.size();i++)
			System.out.println("posMultiSeq: " + posMultiSeq.get(i));
		// // find dNTPs from variation position list by assigned position
		int pos = 32;
		String dntps = "";
		for(int i=0;i<var_pos_list.size();i++) {
			int pos_tmp = var_pos_list.get(i);
			if(pos_tmp == pos) {
				dntps = dntps_list.get(i);
				break;
			}
		}
		System.out.println("dntps: " + dntps);
		// get posMultiSeq from dNTPs seq for full sequence
		List<String> posMultiSeq_dNTPs = seqProcess.getPosMultiSeq_dNTPs(iupac_seq, pos, dntps);
		for(int i=0;i<posMultiSeq_dNTPs.size();i++)
			System.out.println("posMultiSeq: " + posMultiSeq_dNTPs.get(i));
	}*/
}