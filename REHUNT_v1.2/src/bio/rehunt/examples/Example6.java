/*
 * Program name: Example6.java
 * Date: 2016/10/05
 * Update: 2017/10/14
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		REHUNT example 6.
 */

package bio.rehunt.examples;

import java.util.*;

import bio.rehunt.thread.JudgeRFLPBatchThread;

/**
 * REHUNT example 6.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 * Description:
 * 		High throughput analysis.
 * 		The multiple sequences can be analyzed by "JudgeRFLPBatchThread" class.
 * 		The function is useful for high throughput analysis.
 */
public class Example6 {
	public static void main(String args[]) {
		// judge variations in multiple sequences if can be recognized by restriction enzymes
		List<String[]> seqList = new LinkedList<String[]>();
		String[] seq_data_1 = new String[4];
		seq_data_1[0] = "rs137853007";	// rs137853007
		seq_data_1[1] = "ACCGAACATACAGCAAGAAACACTTT";
		seq_data_1[2] = "[C/T]";
		seq_data_1[3] = "GGATTTTCAGGGTAGGTAATGAATA";
		String[] seq_data_2 = new String[4];
		seq_data_2[0] = "rs202217267";	// rs202217267
		seq_data_2[1] = "TGGGACGGCAAGGGGGACTGTAGAT";
		seq_data_2[2] = "[A/G]";
		seq_data_2[3] = "GGTGAAAAGAGCAGTCAGAGGACCA";
		String[] seq_data_3 = new String[4];
		seq_data_3[0] = "rs201930255";	// rs201930255
		seq_data_3[1] = "GCTGGGGCACAGCAGGCCAGTGTGCA";
		seq_data_3[2] = "[C/G]";
		seq_data_3[3] = "GGTGGCAAGTGGCTCCTGACCTGGA";
		seqList.add(seq_data_1);
		seqList.add(seq_data_2);
		seqList.add(seq_data_3);
		JudgeRFLPBatchThread judgeRFLPBatchThread = new JudgeRFLPBatchThread(seqList);
		judgeRFLPBatchThread.start();
		boolean is_wait = true;
		while(judgeRFLPBatchThread.getState() != Thread.State.TERMINATED) {
			if(is_wait) {
				System.out.println("Please waiting...");
				is_wait = false;
			}
		}
		List<Boolean> resultList = judgeRFLPBatchThread.getIsCutResultList();
		for(int i=0;i<resultList.size();i++) {
			boolean isCut = resultList.get(i).booleanValue();
			System.out.println("+ strand for variation " + (i+1) + ": " + isCut);
		}
		List<List<List<String>>> cutEnzymesBatch = judgeRFLPBatchThread.getCutEnzymesBatch();
		for(int i=0;i<cutEnzymesBatch.size();i++) {
			List<List<String>> cutEnzymesList = cutEnzymesBatch.get(i);
			System.out.println("Restriction enzymes information for variation " + (i+1) + ": ");
			System.out.println("Restriction enzymes which can recognize variation1, but cannot recognize variation2: " + cutEnzymesList.get(0).toString());
			System.out.println("Restriction enzymes which can recognize variation2, but cannot recognize variation1: " + cutEnzymesList.get(1).toString());
		}
	}
}