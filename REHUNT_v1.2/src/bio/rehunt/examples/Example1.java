/*
 * Program name: Example1.java
 * Date: 2016/10/05
 * Update: 2017/10/14
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		REHUNT example 1.
 */

package bio.rehunt.examples;

import java.util.*;

import bio.rehunt.seq.Sequence;
import bio.rehunt.rflp.RFLPprocess;
import bio.rehunt.rflp.JudgeRFLP;

/**
 * REHUNT example 1.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 * Description:
 * 		Universal judgement of whether a variation of a sequence can be recognized by restriction enzymes.
 * 		The input variations may include two variations, three variations, or four variations.
 */
public class Example1 {
	public static void main(String args[]) {
		// judge the variation if can be recognized by restriction enzymes
		//String seq = "TTAGCATCAGCATTTGCTGC[A/G/C/T]ATCGCTAACGGTGGATCTAC";
		String seq ="GCGGCACGAGCAGACCCCTGTGTGCCGTCCTGTGGGCGCGGGGCGGCAGGGGAGGCGCACACCTGCTCCTTTGTGCAGCCTCCCCCCTCCCGCAAAGTTAAAGAGCAGGAAAGTCAGGATTCCTCGCTCGGCCCTGCCCTGCCGGCTGCTCCGCGCTCCGCTCCTCCCTGCGAGCGTGTGTGTGTGTCGGGGGTCCCTCCCCTCCTGGCTCTGGGGTCGGGCGCGCACCCCGCCCCGTAGCGCGGCCCCTCCCTGGCGAGCGCAACCCCATCCAGCGGGAGCGCGGAGCCGCGGCCGCGGGGAAGCATTAAGTTTATTCGCCTCAAAGTGACGCAAAAATTCTTCAAGAGCTCTTTGGCGGCGGCTATCTAGAGATCAGACCATGTGAGGGCCCGCGGGTACAAATACGGCCGCGCCGGCGCCCCTCCGCACAGCCAGCGCCGCCGGGTGCCTCGAGGGCGCGAGGCCAGCCCGCCTGCCCAGCCCGGGACCAGCCTCCC[C/T]GCGCAGCCTGGCAGGTGGGTCCGCTTTTCCTCTCCGCCTCGAACCCACGTTTCTTTCCAGACCTTCTTCCCCGCCTCGGGGAGGGGGATAGAACCGCTGCGCCCCACCGCCCTGCGAGGAGGCGAGGAGGTGCATGCGCCCCAGCGGTGGGCGCCGGATCCTGCCCCTGCGCCCTCCACGCTCAGCAAGAGCCAGAGCTGAAGCTGACCGGCCAGAGTGGGAGACGAGGAACGTGGAGTGCTCGAAGTGGGCGGGCGTAGGGGGCTCCTTTGTCTATTGTTGCAGGGGCTTTGCAACCTCTTGAGGACCCGTGAAAGACCCCTTAAAGAGGGTCTTGAGCAAAGTGCTGTCCTGCGCTTTGGGAAGAGTTGCTTGCTTTCGTTTCAACCCATGGTTTATGTTTGATCTTTACTTTGCTGTCATCGCGTGCAGGTGGCTTTATGCAAAGGGAGAGTTCTGGTTGACACAAATGCCCAGACAGCTAGAGAAATCTCTGAGTG";
		RFLPprocess rflpProcess = new RFLPprocess();
		// Sequence
		Sequence sequence = new Sequence();
		String seq_complementary = sequence.complementaryTrans(seq);
		List<List<String>> cutEnzymesList = null;
		List<List<String>> cutEnzymesList_complementary = null;
		// ----- variation1 and variation2 ----- //
		// cutEnzymesList.get(0): (variation1 and variation2), cutEnzymesList.get(1): (variation2 and variation1)
		// ----- variation1, variation2 and variation3 ----- //
		// cutEnzymesList.get(0): (variation1 and variation2), cutEnzymesList.get(1): (variation1 and variation3)
		// cutEnzymesList.get(2): (variation2 and variation1), cutEnzymesList.get(3): (variation2 and variation3)
		// cutEnzymesList.get(4): (variation3 and variation1), cutEnzymesList.get(5): (variation3 and variation2)
		// ----- variation1, variation2, variation3 and variation4 ----- //
		// cutEnzymesList.get(0): (variation1 and variation2), cutEnzymesList.get(1): (variation1 and variation3), cutEnzymesList.get(2): (variation1 and variation4)
		// cutEnzymesList.get(3): (variation2 and variation1), cutEnzymesList.get(4): (variation2 and variation3), cutEnzymesList.get(5): (variation2 and variation4)
		// cutEnzymesList.get(6): (variation3 and variation1), cutEnzymesList.get(7): (variation3 and variation2), cutEnzymesList.get(8): (variation3 and variation4)
		// cutEnzymesList.get(9): (variation4 and variation1), cutEnzymesList.get(10): (variation4 and variation2), cutEnzymesList.get(11): (variation4 and variation3)
		// do RFLP for the sequence
		JudgeRFLP judgeRFLP = new JudgeRFLP(seq);
		judgeRFLP.setIUPACenzyme(true);
		if(judgeRFLP.isVarSeq()) {
			boolean isCut = judgeRFLP.isCanCut_dNTPs();
			System.out.println("+ strand: " + isCut);
			cutEnzymesList = judgeRFLP.getCutEnzymesList();
			System.out.println("cutEnzymesList: " + cutEnzymesList.size());
			// enzymes
			if(cutEnzymesList.size() == 2) {
				System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2: " + cutEnzymesList.get(0).toString());
				System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1: " + cutEnzymesList.get(1).toString());
			}
			else if(cutEnzymesList.size() == 6) {
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
			else if(cutEnzymesList.size() == 12) {
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
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");

		// do RFLP for complementary sequence
		JudgeRFLP judgeRFLP_complementary = new JudgeRFLP(seq_complementary);
		if(judgeRFLP_complementary.isVarSeq()) {
			boolean isCut = judgeRFLP_complementary.isCanCut_dNTPs();
			System.out.println("- strand: " + isCut);
			cutEnzymesList_complementary = judgeRFLP_complementary.getCutEnzymesList();
			System.out.println("cutEnzymesList_complementary: " + cutEnzymesList_complementary.size());
			// enzymes for the complementary sequence
			if(cutEnzymesList_complementary.size() == 2) {
				System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2: " + cutEnzymesList_complementary.get(0).toString());
				System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1: " + cutEnzymesList_complementary.get(1).toString());
			}
			else if(cutEnzymesList_complementary.size() == 6) {
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
			else if(cutEnzymesList_complementary.size() == 12) {
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
		}
		else	// not variation sequence
			System.out.println("The sequence is not a variation sequence.");
	}
}