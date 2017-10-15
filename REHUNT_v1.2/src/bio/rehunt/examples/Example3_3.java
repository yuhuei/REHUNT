/*
 * Program name: Example3_3.java
 * Date: 2016/10/05
 * Update: 2017/10/14
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		REHUNT example 3_3.
 */

package bio.rehunt.examples;

import java.util.*;

import bio.rehunt.seq.Sequence;
import bio.rehunt.rflp.RFLPprocess;
import bio.rehunt.rflp.JudgeRFLP;

/**
 * REHUNT example 3_3.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 * Description:
 * 		Identify restriction enzymes for a sequence and complementary sequence with four variations.
 * 		The sequence "TTAGCATCAGCATTTGCTGC[A/G/C/T]ATCGCTAACGGTGGATCTAC" with four variations A, G, C and T that can be recognized by restriction enzymes.
 * 		Its complementary sequence is also can be recognized by restriction enzymes.
 * 		These restriction enzymes are easy found out by REHUNT.
 */
public class Example3_3 {
	public static void main(String args[]) {
		// judge the variation if can be recognized by restriction enzymes
		String seq = "TTAGCATCAGCATTTGCTGC[A/G/C/T]ATCGCTAACGGTGGATCTAC";
		RFLPprocess rflpProcess = new RFLPprocess();
		// Sequence
		Sequence sequence = new Sequence();
		String seq_complementary = sequence.complementaryTrans(seq);
		List<List<String>> cutEnzymesList = null;
		List<List<String>> cutEnzymesList_complementary = null;
		// ----- variation1, variation2, variation3 and variation4 ----- //
		// cutEnzymesList.get(0): (variation1 and variation2), cutEnzymesList.get(1): (variation1 and variation3), cutEnzymesList.get(2): (variation1 and variation4)
		// cutEnzymesList.get(3): (variation2 and variation1), cutEnzymesList.get(4): (variation2 and variation3), cutEnzymesList.get(5): (variation2 and variation4)
		// cutEnzymesList.get(6): (variation3 and variation1), cutEnzymesList.get(7): (variation3 and variation2), cutEnzymesList.get(8): (variation3 and variation4)
		// cutEnzymesList.get(9): (variation4 and variation1), cutEnzymesList.get(10): (variation4 and variation2), cutEnzymesList.get(11): (variation4 and variation3)
		// do RFLP for sequence with four variations
		JudgeRFLP judgeRFLP = new JudgeRFLP(seq);
		if(judgeRFLP.isVarSeq()) {
			boolean isCut = judgeRFLP.isCanCut_dNTPs();
			System.out.println("+ strand: " + isCut);
			cutEnzymesList = judgeRFLP.getCutEnzymesList();
			// enzymes
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2: " + cutEnzymesList.get(0).toString());
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation3: " + cutEnzymesList.get(1).toString());
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation4: " + cutEnzymesList.get(2).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1: " + cutEnzymesList.get(3).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation3: " + cutEnzymesList.get(4).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation4: " + cutEnzymesList.get(5).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation1: " + cutEnzymesList.get(6).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation2: " + cutEnzymesList.get(7).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation4: " + cutEnzymesList.get(8).toString());
			System.out.println("Restriction enzymes which can recognize variation4, but cannot recognize variation1: " + cutEnzymesList.get(9).toString());
			System.out.println("Restriction enzymes which can recognize variation4, but cannot recognize variation2: " + cutEnzymesList.get(10).toString());
			System.out.println("Restriction enzymes which can recognize variation4, but cannot recognize variation3: " + cutEnzymesList.get(11).toString());
			List<String> cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymesList.get(0), cutEnzymesList.get(1));
			cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymes_seq1, cutEnzymesList.get(2));
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2, variation3 and variation4: " + cutEnzymes_seq1);
			List<String> cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymesList.get(3), cutEnzymesList.get(4));
			cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymes_seq2, cutEnzymesList.get(5));
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1, variation3 and variation4: " + cutEnzymes_seq2);
			List<String> cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymesList.get(6), cutEnzymesList.get(7));
			cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymes_seq3, cutEnzymesList.get(8));
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation1, variation2 and variation4: " + cutEnzymes_seq3);
			List<String> cutEnzymes_seq4 = rflpProcess.getSame(cutEnzymesList.get(9), cutEnzymesList.get(10));
			cutEnzymes_seq4 = rflpProcess.getSame(cutEnzymes_seq4, cutEnzymesList.get(11));
			System.out.println("Restriction enzymes which can recognize variation4, but cannot recognize variation1, variation2 and variation3: " + cutEnzymes_seq4);
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");

		// do RFLP for complementary sequence with four variations
		JudgeRFLP judgeRFLP_complementary = new JudgeRFLP(seq_complementary);
		if(judgeRFLP_complementary.isVarSeq()) {
			boolean isCut = judgeRFLP_complementary.isCanCut_dNTPs();
			System.out.println("- strand: " + isCut);
			cutEnzymesList_complementary = judgeRFLP_complementary.getCutEnzymesList();
			// enzymes for the complementary sequence
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2: " + cutEnzymesList_complementary.get(0).toString());
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation3: " + cutEnzymesList_complementary.get(1).toString());
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation4: " + cutEnzymesList_complementary.get(2).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1: " + cutEnzymesList_complementary.get(3).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation3: " + cutEnzymesList_complementary.get(4).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation4: " + cutEnzymesList_complementary.get(5).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation1: " + cutEnzymesList_complementary.get(6).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation2: " + cutEnzymesList_complementary.get(7).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation4: " + cutEnzymesList_complementary.get(8).toString());
			System.out.println("Restriction enzymes which can recognize variation4, but cannot recognize variation1: " + cutEnzymesList_complementary.get(9).toString());
			System.out.println("Restriction enzymes which can recognize variation4, but cannot recognize variation2: " + cutEnzymesList_complementary.get(10).toString());
			System.out.println("Restriction enzymes which can recognize variation4, but cannot recognize variation3: " + cutEnzymesList_complementary.get(11).toString());
			List<String> cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymesList_complementary.get(0), cutEnzymesList_complementary.get(1));
			cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymes_seq1, cutEnzymesList_complementary.get(2));
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2, variation3 and variation4: " + cutEnzymes_seq1);
			List<String> cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymesList_complementary.get(3), cutEnzymesList_complementary.get(4));
			cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymes_seq2, cutEnzymesList_complementary.get(5));
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1, variation3 and variation4: " + cutEnzymes_seq2);
			List<String> cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymesList_complementary.get(6), cutEnzymesList_complementary.get(7));
			cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymes_seq3, cutEnzymesList_complementary.get(8));
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation1, variation2 and variation4: " + cutEnzymes_seq3);
			List<String> cutEnzymes_seq4 = rflpProcess.getSame(cutEnzymesList_complementary.get(9), cutEnzymesList_complementary.get(10));
			cutEnzymes_seq4 = rflpProcess.getSame(cutEnzymes_seq4, cutEnzymesList_complementary.get(11));
			System.out.println("Restriction enzymes which can recognize variation4, but cannot recognize variation1, variation2 and variation3: " + cutEnzymes_seq4);
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");
	}
}