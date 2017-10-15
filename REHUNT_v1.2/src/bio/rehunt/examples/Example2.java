/*
 * Program name: Example2.java
 * Date: 2016/10/05
 * Update: 2017/10/14
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		REHUNT example 2.
 */

package bio.rehunt.examples;

import java.util.*;

import bio.rehunt.rflp.JudgeRFLP;

/**
 * REHUNT example 2.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 * Description:
 * 		Search for specific restriction enzymes.
 * 		For example, the sequence "ACGG[A/C]TTTTTT" can be recognized by restriction enzyme TspGWI (ACGGA) for variation A.
 * 		The sequence "ACGG[A/C]TTTTTTACGGATTT" can be recognized by restriction enzyme TspGWI (ACGGA) for variation A, but it will be excluded because of the repeat of the sequence "ACGGA".
 * 		REHUNT identifies specific and available restriction enzymes, thus, the reappearing restriction enzyme TspGWI (ACGGA) will be excluded.
 */
public class Example2 {
	public static void main(String args[]) {
		// judge if the variation can be recognized by restriction enzymes
		String seq1 = "ACGG[A/C]TTTTTT";	// restriction enzyme TspGWI (ACGGA) can recognize variation A
		String seq2 = "ACGG[A/C]TTTTTTACGGATTT";	// restriction enzyme TspGWI (ACGGA) can recognize variation A, but it is excluded
		List<List<String>> cutEnzymesList = null;
		// ----- variation1 and variation2 ----- //
		// cutEnzymesList.get(0): (restriction enzymes which can recognize variation1, but cannot recognize variation2)
		// do RFLP for seq1
		JudgeRFLP judgeRFLP = new JudgeRFLP(seq1);
		if(judgeRFLP.isVarSeq()) {
			// find specific restriction enzymes for recognizing variations
			boolean isCut = judgeRFLP.isCanCut_dNTPs();
			System.out.println("+ strand: " + isCut);
			if(isCut) {
				// enzymes
				cutEnzymesList = judgeRFLP.getCutEnzymesList();
				System.out.println("seq1 - Restriction enzymes which can recognize variation A, but cannot recognize variation C: " + cutEnzymesList.get(0).toString());
			}
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");
		
		// do RFLP for seq2
		judgeRFLP.setSeq(seq2);
		if(judgeRFLP.isVarSeq()) {
			// find specific restriction enzymes for recognizing variations
			boolean isCut = judgeRFLP.isCanCut_dNTPs();
			System.out.println("+ strand: " + isCut);
			if(isCut) {
				// enzymes
				cutEnzymesList = judgeRFLP.getCutEnzymesList();
				System.out.println("seq2 - Restriction enzymes which can recognize variation A, but cannot recognize variation C: " + cutEnzymesList.get(0).toString());
			}
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");
	}
}