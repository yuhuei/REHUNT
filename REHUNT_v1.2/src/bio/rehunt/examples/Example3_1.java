/*
 * Program name: Example3_1.java
 * Date: 2016/10/05
 * Update: 2017/10/14
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		REHUNT example 3_1.
 */

package bio.rehunt.examples;

import java.util.*;

import bio.rehunt.seq.Sequence;
import bio.rehunt.rflp.JudgeRFLP;

/**
 * REHUNT example 3_1.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 * Description:
 * 		Identify restriction enzymes for a sequence and complementary sequence with two variations.
 * 		The sequence "TTAGCATCAGCATTTGCTGC[A/G]ATCGCTAACGGTGGATCTAC" with two variations A and G that can be recognized by restriction enzymes.
 * 		Its complementary sequence is also can be recognized by restriction enzymes.
 * 		These restriction enzymes are easy found out by REHUNT.
 */
public class Example3_1 {
	public static void main(String args[]) {
		// judge the variation if can be recognized by restriction enzymes
		String seq = "TTAGCATCAGCATTTGCTGC[A/G]ATCGCTAACGGTGGATCTAC";
		// Sequence
		Sequence sequence = new Sequence();
		String seq_complementary = sequence.complementaryTrans(seq);
		List<List<String>> cutEnzymesList = null;
		List<List<String>> cutEnzymesList_complementary = null;
		// ----- variation1 and variation2 ----- //
		// cutEnzymesList.get(0): (variation1 and variation2), cutEnzymesList.get(1): (variation2 and variation1)
		// do RFLP for sequence with two variations
		JudgeRFLP judgeRFLP = new JudgeRFLP(seq);
		if(judgeRFLP.isVarSeq()) {
			boolean isCut = judgeRFLP.isCanCut_dNTPs();
			System.out.println("+ strand: " + isCut);
			cutEnzymesList = judgeRFLP.getCutEnzymesList();
			// enzymes
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2: " + cutEnzymesList.get(0).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1: " + cutEnzymesList.get(1).toString());
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");

		// do RFLP for complementary sequence with two variations
		JudgeRFLP judgeRFLP_complementary = new JudgeRFLP(seq_complementary);
		if(judgeRFLP_complementary.isVarSeq()) {
			boolean isCut = judgeRFLP_complementary.isCanCut_dNTPs();
			System.out.println("- strand: " + isCut);
			cutEnzymesList_complementary = judgeRFLP_complementary.getCutEnzymesList();
			// enzymes for the complementary sequence
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2: " + cutEnzymesList_complementary.get(0).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1: " + cutEnzymesList_complementary.get(1).toString());
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");
	}
}