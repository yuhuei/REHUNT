/*
 * Program name: JudgeRFLP.java
 * Date: 2016/09/30
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		Judge if variation can be recognized by restriction enzymes.
 */

package bio.rehunt.rflp;

import java.util.*;

import bio.rehunt.seq.Sequence;
import bio.rehunt.seq.SeqProcess;

/**
 * Judge if variation can be recognized by restriction enzymes.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0_
 */
public class JudgeRFLP {
	private int var_pos = -1;	// default var_pos, -1 is the first variation position
	private int enzymeLen_min = 4;
	private int enzymeLen_max = 10;
	private boolean isIUPACenzyme = false;	// judge whether get IUPAC enzymes or not
	private boolean is_var_seq = false;
	private String iupac_seq = null;
	private List<String> dntps_list = null;
	private List<Integer> var_pos_list = null;
	private List<List<String>> enzymeNameList = null;
	private List<List<String>> cutEnzymesList = null;

	/**
	 * Constructor for JudgeRFLP.
	 */
	public JudgeRFLP() {}
	
	/**
	 * Constructor for initialize JudgeRFLP.
	 * @param seq Sequence.
	 */
	public JudgeRFLP(String seq) {
		setSeq(seq);
	}
	
	/**
	 * Judge if the sequence has variations.
	 * @return If variation sequence then return true, else return false.
	 */
	public boolean isVarSeq() {
		return is_var_seq;
	}
	
	/**
	 * Set sequence.
	 * @param seq Sequence.
	 */
	public void setSeq(String seq) {
		Sequence sequence = new Sequence();
		// judge seq if variation sequence and get related information
		if(is_var_seq=sequence.isVarSeq(seq)) {
			boolean is_success = sequence.makeVarSeq(seq);
			if(is_success) {
				iupac_seq = sequence.getIUPACSeq();
				dntps_list = sequence.getdNTPsList();
				var_pos_list = sequence.getVarPosList();
			}
		}
		this.enzymeNameList = new LinkedList<List<String>>();
		this.cutEnzymesList = new LinkedList<List<String>>();
	}

	/**
	 * Set variation position in sequence for RFLP analysis.
	 * @param var_pos Variation position.
	 */
	public void setVarPos(int var_pos) {
		this.var_pos = var_pos;
	}

	/**
	 * Set if to find the IUPAC enzymes.
	 * @param isIUPACenzyme True is to find the IUPAC enzymes, and false is not.
	 */
	public void setIUPACenzyme(boolean isIUPACenzyme) {
		this.isIUPACenzyme = isIUPACenzyme;
	}

	/**
	 * Set minimum length of enzyme sequence for search.
	 * @param enzyme_len Restriction enzyme length.
	 */
	public void setEnzymeLenMin(int enzyme_len) {
		enzymeLen_min = enzyme_len;
	}

	/**
	 * Set maximum length of enzyme sequence for search.
	 * @param enzyme_len Restriction enzyme length.
	 */
	public void setEnzymeLenMax(int enzyme_len) {
		enzymeLen_max = enzyme_len;
	}

	/**
	 * Judge if sequence of IUPAC format can be recognized by restriction enzymes.
	 * @return If return true then sequence can be recognized by restriction enzymes, else return false.
	 */
	public boolean isCanCut() {
		// clear
		enzymeNameList.clear();
		cutEnzymesList.clear();

		boolean isCut = false;
		// SeqProcess
		SeqProcess seqProcess = new SeqProcess();
		if(var_pos == -1)	// if var_pos=-1, set var_pos is first IUPAC position
			var_pos = seqProcess.getFirstIUPACPos(iupac_seq);
		String seq = iupac_seq;
		int seq_var_pos = var_pos;
		// get posMultiSeq from IUPAC seq for full sequence
		List<String> posMultiSeq = seqProcess.getPosMultiSeq(seq, seq_var_pos);
		/*for(int i=0;i<posMultiSeq.size();i++)
			System.out.println("posMultiSeq: " + posMultiSeq.get(i));*/
		// RFLP
		RFLP rflp = new RFLP();
		rflp.setIUPACenzyme(isIUPACenzyme);
		rflp.setEnzymeLenMin(enzymeLen_min);
		rflp.setEnzymeLenMax(enzymeLen_max);
		for(int i=0;i<posMultiSeq.size();i++) {
			String singleSeq = posMultiSeq.get(i);
			rflp.setSeq(singleSeq);
			if(rflp.getEnzymes()) {
				try {
					List<String> enzymeName_list = rflp.getEnzymeNameList();
					enzymeNameList.add(enzymeName_list);
					/*for(int j=0;j<enzymeName_list.size();j++)
						System.out.println(enzymeName_list.get(j));*/
				}
				catch(Exception e) {
					System.out.println("Exception: " + e);
				}
			}
			else
				System.out.println("Error: Not to get Enzymes.");
		}
		// judge seq if can be recognized by restriction enzymes
		RFLPprocess rflpProcess = new RFLPprocess();
		List<String> enzymeNameList0 = enzymeNameList.get(0);
		for(int i=1;i<enzymeNameList.size();i++) {
			List<String> enzymeName_list = enzymeNameList.get(i);
			List<String> cutEnzymes = rflpProcess.getDiffOrder(enzymeNameList0, enzymeName_list);
			if(cutEnzymes.size() != 0) {
				isCut = true;
				//System.out.println("This sequence can be recognized.");
				break;
			}
			else {
				//System.out.println("This sequence can't be recognized.");
				isCut = false;
			}
		}
		// find cutEnzymes cut which singleSeq
		for(int i=0;i<enzymeNameList.size();i++) {
			List<String> enzymeName_listA = enzymeNameList.get(i);
			for(int j=0;j<enzymeNameList.size();j++) {
				List<String> enzymeName_listB = enzymeNameList.get(j);
				if(i != j) {
					List<String> cutEnzymes = rflpProcess.getDiff(enzymeName_listA, enzymeName_listB);
					cutEnzymesList.add(cutEnzymes);
				}
			}
		}
		return isCut;
	}

	/**
	 * Judge if sequence of dNTPs format can be recognized by restriction enzymes.
	 * @return If return true then sequence can be recognized by restriction enzymes, else return false.
	 */
	public boolean isCanCut_dNTPs() {
		// clear
		enzymeNameList.clear();
		cutEnzymesList.clear();

		boolean isCut = false;
		// SeqProcess
		SeqProcess seqProcess = new SeqProcess();
		if(var_pos == -1)	// if var_pos=-1, set var_pos is first IUPAC position
			var_pos = seqProcess.getFirstIUPACPos(iupac_seq);
		String seq = iupac_seq;
		int seq_var_pos = var_pos;
		// find dntps where position on var_pos
		String dntps = "";
		for(int i=0;i<var_pos_list.size();i++) {
			int pos = var_pos_list.get(i);
			if(var_pos == pos) {
				dntps = dntps_list.get(i);
				break;
			}
		}
		// get posMultiSeq from dNTPs seq for full sequence
		List<String> posMultiSeq_dNTPs = seqProcess.getPosMultiSeq_dNTPs(seq, seq_var_pos, dntps);
		/*for(int i=0;i<posMultiSeq_dNTPs.size();i++)
			System.out.println("posMultiSeq: " + posMultiSeq_dNTPs.get(i));*/
		// RFLP
		RFLP rflp = new RFLP();
		rflp.setIUPACenzyme(isIUPACenzyme);
		rflp.setEnzymeLenMin(enzymeLen_min);
		rflp.setEnzymeLenMax(enzymeLen_max);
		for(int i=0;i<posMultiSeq_dNTPs.size();i++) {
			String singleSeq = posMultiSeq_dNTPs.get(i);
			rflp.setSeq(singleSeq);
			if(rflp.getEnzymes()) {
				try {
					List<String> enzymeName_list = rflp.getEnzymeNameList();
					enzymeNameList.add(enzymeName_list);
					/*for(int j=0;j<enzymeName_list.size();j++)
						System.out.println(enzymeName_list.get(j));*/
				}
				catch(Exception e) {
					System.out.println("Exception: " + e);
				}
			}
			else
				System.out.println("Error: Fail to get Enzymes.");
		}
		// judge seq if can be cut
		RFLPprocess rflpProcess = new RFLPprocess();
		List<String> enzymeNameList0 = enzymeNameList.get(0);
		for(int i=1;i<enzymeNameList.size();i++) {
			List<String> enzymeName_list = enzymeNameList.get(i);
			List<String> cutEnzymes = rflpProcess.getDiffOrder(enzymeNameList0, enzymeName_list);
			if(cutEnzymes.size() != 0) {
				isCut = true;
				//System.out.println("This sequence can be recognized.");
				break;
			}
			else {
				//System.out.println("This sequence can't be recognized.");
				isCut = false;
			}
		}
		// find cutEnzymes cut which singleSeq
		for(int i=0;i<enzymeNameList.size();i++) {
			List<String> enzymeName_listA = enzymeNameList.get(i);
			for(int j=0;j<enzymeNameList.size();j++) {
				List<String> enzymeName_listB = enzymeNameList.get(j);
				if(i != j) {
					List<String> cutEnzymes = rflpProcess.getDiff(enzymeName_listA, enzymeName_listB);
					cutEnzymesList.add(cutEnzymes);
				}
			}
		}
		return isCut;
	}

	/**
	 * Get restriction enzyme name list that can recognize variation.
	 * @return Restriction enzyme name list that can recognize variation.
	 */
	public List<List<String>> getCutEnzymesList() {
		// ----- seq1 and seq2 ----- //
		// cutEnzymesList.get(0): (seq1 and seq2), cutEnzymesList.get(1): (seq2 and seq1)
		// ----- seq1, seq2 and seq3 ----- //
		// cutEnzymesList.get(0): (seq1 and seq2), cutEnzymesList.get(1): (seq1 and seq3)
		// cutEnzymesList.get(2): (seq2 and seq1), cutEnzymesList.get(3): (seq2 and seq3)
		// cutEnzymesList.get(4): (seq3 and seq1), cutEnzymesList.get(5): (seq3 and seq2)
		// ----- seq1, seq2, seq3 and seq4 ----- //
		// cutEnzymesList.get(0): (seq1 and seq2), cutEnzymesList.get(1): (seq1 and seq3), cutEnzymesList.get(2): (seq1 and seq4)
		// cutEnzymesList.get(3): (seq2 and seq1), cutEnzymesList.get(4): (seq2 and seq3), cutEnzymesList.get(5): (seq2 and seq4)
		// cutEnzymesList.get(6): (seq3 and seq1), cutEnzymesList.get(7): (seq3 and seq2), cutEnzymesList.get(8): (seq3 and seq4)
		// cutEnzymesList.get(9): (seq4 and seq1), cutEnzymesList.get(10): (seq4 and seq2), cutEnzymesList.get(11): (seq4 and seq3)
		return cutEnzymesList;
	}

	/**
	 * Get restriction enzyme name list that exist in the sequence.
	 * @return Restriction enzyme name list that exist in the sequence.
	 */
	// get enzymeNameList
	public List<List<String>> getEnzymeNameList() {
		return enzymeNameList;
	}

/*	public static void main(String args[]) {
		// judge sequence have restriction enzyme
		String seq = "TATTCAAGTGCACGAGACCAATGAC[A/C/G/T]GGACCTCTGGTGAGGCCCTGGTGAG";
		RFLPprocess rflpProcess = new RFLPprocess();
		// Sequence
		Sequence sequence = new Sequence();
		String seq_complementary = sequence.complementaryTrans(seq);
		List<List<String>> cutEnzymesList = null;
		List<List<String>> cutEnzymesList_complementary = null;
		// ----- seq1 and seq2 ----- //
		// cutEnzymesList.get(0): (seq1 and seq2), cutEnzymesList.get(1): (seq2 and seq1)
		// ----- seq1, seq2 and seq3 ----- //
		// cutEnzymesList.get(0): (seq1 and seq2), cutEnzymesList.get(1): (seq1 and seq3)
		// cutEnzymesList.get(2): (seq2 and seq1), cutEnzymesList.get(3): (seq2 and seq3)
		// cutEnzymesList.get(4): (seq3 and seq1), cutEnzymesList.get(5): (seq3 and seq2)
		// ----- seq1, seq2, seq3 and seq4 ----- //
		// cutEnzymesList.get(0): (seq1 and seq2), cutEnzymesList.get(1): (seq1 and seq3), cutEnzymesList.get(2): (seq1 and seq4)
		// cutEnzymesList.get(3): (seq2 and seq1), cutEnzymesList.get(4): (seq2 and seq3), cutEnzymesList.get(5): (seq2 and seq4)
		// cutEnzymesList.get(6): (seq3 and seq1), cutEnzymesList.get(7): (seq3 and seq2), cutEnzymesList.get(8): (seq3 and seq4)
		// cutEnzymesList.get(9): (seq4 and seq1), cutEnzymesList.get(10): (seq4 and seq2), cutEnzymesList.get(11): (seq4 and seq3)
		// do seq RFLP for variations
		JudgeRFLP judgeRFLP = new JudgeRFLP(seq);
		if(judgeRFLP.isVarSeq()) {
			boolean isCut = judgeRFLP.isCanCut_dNTPs();
			System.out.println("+ strand: " + isCut);
			cutEnzymesList = judgeRFLP.getCutEnzymesList();
			System.out.println("cutEnzymesList: " + cutEnzymesList.size());
			// enzymes
			if(cutEnzymesList.size() == 2) {
				System.out.println("cut seq1 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList.get(0).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList.get(1).toString());
			}
			else if(cutEnzymesList.size() == 6) {
				System.out.println("cut seq1 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList.get(0).toString());
				System.out.println("cut seq1 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList.get(1).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList.get(2).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList.get(3).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList.get(4).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList.get(5).toString());
				List<String> cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymesList.get(0), cutEnzymesList.get(1));
				System.out.println("cut seq1 Enzymes, but can't cut seq2 and seq3 Enzymes: " + cutEnzymes_seq1);
				List<String> cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymesList.get(2), cutEnzymesList.get(3));
				System.out.println("cut seq2 Enzymes, but can't cut seq1 and seq3 Enzymes: " + cutEnzymes_seq2);
				List<String> cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymesList.get(4), cutEnzymesList.get(5));
				System.out.println("cut seq3 Enzymes, but can't cut seq1 and seq2 Enzymes: " + cutEnzymes_seq3);
			}
			else if(cutEnzymesList.size() == 12) {
				System.out.println("cut seq1 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList.get(0).toString());
				System.out.println("cut seq1 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList.get(1).toString());
				System.out.println("cut seq1 Enzymes, but can't cut seq4 Enzymes: " + cutEnzymesList.get(2).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList.get(3).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList.get(4).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq4 Enzymes: " + cutEnzymesList.get(5).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList.get(6).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList.get(7).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq4 Enzymes: " + cutEnzymesList.get(8).toString());
				System.out.println("cut seq4 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList.get(9).toString());
				System.out.println("cut seq4 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList.get(10).toString());
				System.out.println("cut seq4 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList.get(11).toString());
				List<String> cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymesList.get(0), cutEnzymesList.get(1));
				cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymes_seq1, cutEnzymesList.get(2));
				System.out.println("cut seq1 Enzymes, but can't cut seq2, seq3 and seq4 Enzymes: " + cutEnzymes_seq1);
				List<String> cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymesList.get(3), cutEnzymesList.get(4));
				cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymes_seq2, cutEnzymesList.get(5));
				System.out.println("cut seq2 Enzymes, but can't cut seq1, seq3 and seq4 Enzymes: " + cutEnzymes_seq2);
				List<String> cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymesList.get(6), cutEnzymesList.get(7));
				cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymes_seq3, cutEnzymesList.get(8));
				System.out.println("cut seq3 Enzymes, but can't cut seq1, seq2 and seq4 Enzymes: " + cutEnzymes_seq3);
				List<String> cutEnzymes_seq4 = rflpProcess.getSame(cutEnzymesList.get(9), cutEnzymesList.get(10));
				cutEnzymes_seq4 = rflpProcess.getSame(cutEnzymes_seq4, cutEnzymesList.get(11));
				System.out.println("cut seq4 Enzymes, but can't cut seq1, seq2 and seq3 Enzymes: " + cutEnzymes_seq4);
			}
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");

		// do seq_complementary RFLP
		JudgeRFLP judgeRFLP_complementary = new JudgeRFLP(seq_complementary);
		if(judgeRFLP_complementary.isVarSeq()) {
			boolean isCut = judgeRFLP_complementary.isCanCut_dNTPs();
			System.out.println("- strand: " + isCut);
			cutEnzymesList_complementary = judgeRFLP_complementary.getCutEnzymesList();
			System.out.println("cutEnzymesList_complementary: " + cutEnzymesList_complementary.size());
			// enzymes_complementary
			if(cutEnzymesList_complementary.size() == 2) {
				System.out.println("cut seq1 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList_complementary.get(0).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList_complementary.get(1).toString());
			}
			else if(cutEnzymesList_complementary.size() == 6) {
				System.out.println("cut seq1 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList_complementary.get(0).toString());
				System.out.println("cut seq1 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList_complementary.get(1).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList_complementary.get(2).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList_complementary.get(3).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList_complementary.get(4).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList_complementary.get(5).toString());
				List<String> cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymesList_complementary.get(0), cutEnzymesList_complementary.get(1));
				System.out.println("cut seq1 Enzymes, but can't cut seq2 and seq3 Enzymes: " + cutEnzymes_seq1);
				List<String> cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymesList_complementary.get(2), cutEnzymesList_complementary.get(3));
				System.out.println("cut seq2 Enzymes, but can't cut seq1 and seq3 Enzymes: " + cutEnzymes_seq2);
				List<String> cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymesList_complementary.get(4), cutEnzymesList_complementary.get(5));
				System.out.println("cut seq3 Enzymes, but can't cut seq1 and seq2 Enzymes: " + cutEnzymes_seq3);
			}
			else if(cutEnzymesList_complementary.size() == 12) {
				System.out.println("cut seq1 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList_complementary.get(0).toString());
				System.out.println("cut seq1 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList_complementary.get(1).toString());
				System.out.println("cut seq1 Enzymes, but can't cut seq4 Enzymes: " + cutEnzymesList_complementary.get(2).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList_complementary.get(3).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList_complementary.get(4).toString());
				System.out.println("cut seq2 Enzymes, but can't cut seq4 Enzymes: " + cutEnzymesList_complementary.get(5).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList_complementary.get(6).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList_complementary.get(7).toString());
				System.out.println("cut seq3 Enzymes, but can't cut seq4 Enzymes: " + cutEnzymesList_complementary.get(8).toString());
				System.out.println("cut seq4 Enzymes, but can't cut seq1 Enzymes: " + cutEnzymesList_complementary.get(9).toString());
				System.out.println("cut seq4 Enzymes, but can't cut seq2 Enzymes: " + cutEnzymesList_complementary.get(10).toString());
				System.out.println("cut seq4 Enzymes, but can't cut seq3 Enzymes: " + cutEnzymesList_complementary.get(11).toString());
				List<String> cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymesList_complementary.get(0), cutEnzymesList_complementary.get(1));
				cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymes_seq1, cutEnzymesList_complementary.get(2));
				System.out.println("cut seq1 Enzymes, but can't cut seq2, seq3 and seq4 Enzymes: " + cutEnzymes_seq1);
				List<String> cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymesList_complementary.get(3), cutEnzymesList_complementary.get(4));
				cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymes_seq2, cutEnzymesList_complementary.get(5));
				System.out.println("cut seq2 Enzymes, but can't cut seq1, seq3 and seq4 Enzymes: " + cutEnzymes_seq2);
				List<String> cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymesList_complementary.get(6), cutEnzymesList_complementary.get(7));
				cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymes_seq3, cutEnzymesList_complementary.get(8));
				System.out.println("cut seq3 Enzymes, but can't cut seq1, seq2 and seq4 Enzymes: " + cutEnzymes_seq3);
				List<String> cutEnzymes_seq4 = rflpProcess.getSame(cutEnzymesList_complementary.get(9), cutEnzymesList_complementary.get(10));
				cutEnzymes_seq4 = rflpProcess.getSame(cutEnzymes_seq4, cutEnzymesList_complementary.get(11));
				System.out.println("cut seq4 Enzymes, but can't cut seq1, seq2 and seq3 Enzymes: " + cutEnzymes_seq4);
			}
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");
	}*/
}