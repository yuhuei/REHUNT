/*
 * Program name: Example3_2.java
 * Date: 2016/10/05
 * Update: 2017/10/14
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		REHUNT example 3_2.
 */

package bio.rehunt.examples;

import java.util.*;

import bio.rehunt.seq.Sequence;
import bio.rehunt.rflp.RFLPprocess;
import bio.rehunt.rflp.JudgeRFLP;

/**
 * REHUNT example 3_2.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 * Description:
 * 		Identify restriction enzymes for a sequence and complementary sequence with three variations.
 * 		The sequence "TTAGCATCAGCATTTGCTGC[A/G/C]ATCGCTAACGGTGGATCTAC" with three variations A, G and C that can be recognized by restriction enzymes.
 * 		Its complementary sequence is also can be recognized by restriction enzymes.
 * 		These restriction enzymes are easy found out by REHUNT.
 */
public class Example3_2 {
	public static void main(String args[]) {
		// judge the variation if can be recognized by restriction enzymes
		String seq = "TTAGCATCAGCATTTGCTGC[A/G/C]ATCGCTAACGGTGGATCTAC";
		RFLPprocess rflpProcess = new RFLPprocess();
		// Sequence
		Sequence sequence = new Sequence();
		String seq_complementary = sequence.complementaryTrans(seq);
		List<List<String>> cutEnzymesList = null;
		List<List<String>> cutEnzymesList_complementary = null;
		// ----- variation1, svariation2 and variation3 ----- //
		// cutEnzymesList.get(0): (variation1 and variation2), cutEnzymesList.get(1): (variation1 and variation3)
		// cutEnzymesList.get(2): (variation2 and variation1), cutEnzymesList.get(3): (variation2 and variation3)
		// cutEnzymesList.get(4): (variation3 and variation1), cutEnzymesList.get(5): (variation3 and variation2)
		// do RFLP for sequence with three variations
		JudgeRFLP judgeRFLP = new JudgeRFLP(seq);
		if(judgeRFLP.isVarSeq()) {
			boolean isCut = judgeRFLP.isCanCut_dNTPs();
			System.out.println("+ strand: " + isCut);
			cutEnzymesList = judgeRFLP.getCutEnzymesList();
			// enzymes
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2: " + cutEnzymesList.get(0).toString());
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation3: " + cutEnzymesList.get(1).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1: " + cutEnzymesList.get(2).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation3: " + cutEnzymesList.get(3).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation1: " + cutEnzymesList.get(4).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation2: " + cutEnzymesList.get(5).toString());
			List<String> cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymesList.get(0), cutEnzymesList.get(1));
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2 and variation3: " + cutEnzymes_seq1);
			List<String> cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymesList.get(2), cutEnzymesList.get(3));
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1 and variation3: " + cutEnzymes_seq2);
			List<String> cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymesList.get(4), cutEnzymesList.get(5));
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation1 and variation2: " + cutEnzymes_seq3);
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");

		// do RFLP for complementary sequence with three variations
		JudgeRFLP judgeRFLP_complementary = new JudgeRFLP(seq_complementary);
		if(judgeRFLP_complementary.isVarSeq()) {
			boolean isCut = judgeRFLP_complementary.isCanCut_dNTPs();
			System.out.println("- strand: " + isCut);
			cutEnzymesList_complementary = judgeRFLP_complementary.getCutEnzymesList();
			// enzymes for the complementary sequence
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2: " + cutEnzymesList_complementary.get(0).toString());
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation3: " + cutEnzymesList_complementary.get(1).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1: " + cutEnzymesList_complementary.get(2).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation3: " + cutEnzymesList_complementary.get(3).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation1: " + cutEnzymesList_complementary.get(4).toString());
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation2: " + cutEnzymesList_complementary.get(5).toString());
			List<String> cutEnzymes_seq1 = rflpProcess.getSame(cutEnzymesList_complementary.get(0), cutEnzymesList_complementary.get(1));
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2 and variation3: " + cutEnzymes_seq1);
			List<String> cutEnzymes_seq2 = rflpProcess.getSame(cutEnzymesList_complementary.get(2), cutEnzymesList_complementary.get(3));
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1 and variation3: " + cutEnzymes_seq2);
			List<String> cutEnzymes_seq3 = rflpProcess.getSame(cutEnzymesList_complementary.get(4), cutEnzymesList_complementary.get(5));
			System.out.println("Restriction enzymes which can recognize variation3, but cannot recognize variation1 and variation2: " + cutEnzymes_seq3);
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");
	}
}