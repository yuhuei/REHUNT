/*
 * Program name: JudgeRFLPBatchThread.java
 * Date: 2016/09/30
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 * 		Judge RFLP batch thread.
 */

package bio.rehunt.thread;

import java.util.*;

import bio.rehunt.seq.Sequence;
import bio.rehunt.rflp.JudgeRFLP;

/**
 * Judge RFLP batch thread.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 */
public class JudgeRFLPBatchThread extends Thread {
	private List<String[]> seqList = null;
	private List<Boolean> resultList = null;	// the result for judge the variation if can be recognized by restriction enzymes.
	private List<Boolean> resultList_complementary = null;	// the complementary result for judge the variation if can be recognized by restriction enzymes.
	private List<List<List<String>>> cutEnzymesBatch = null;
	private List<List<List<String>>> cutEnzymesBatch_complementary = null;
	private long totalTime = 0l;

	/**
	 * Constructor for initialize JudgeRFLPBatchThread.
	 * @param seqList Sequence list.
	 */
	public JudgeRFLPBatchThread(List<String[]> seqList) {
		this.seqList = seqList;
		resultList = new LinkedList<Boolean>();
		resultList_complementary = new LinkedList<Boolean>();
		cutEnzymesBatch = new LinkedList<List<List<String>>>();
		cutEnzymesBatch_complementary = new LinkedList<List<List<String>>>();
	}

	public void run() {
		// Get current time
		long start = System.currentTimeMillis();
		//----------------------------------------------------------------------
		for(int i=0;i<seqList.size();i++) {
			// get seq_data
			String[] seq_data = seqList.get(i);
			//String seq_id = seq_data[0];
			String sequence5 = seq_data[1];
			String dntps = seq_data[2];
			String sequence3 = seq_data[3];
			
			// get target variation flanking sequence 20 bps
			int flank_len = 20;
			if(sequence5.length() > flank_len)
				sequence5 = sequence5.substring(sequence5.length() - flank_len, sequence5.length());
			if(sequence3.length() > flank_len)
				sequence3 = sequence3.substring(0, flank_len);
			String seq = sequence5 + dntps + sequence3;	// dNTP sequence

			// Sequence
			Sequence sequence = new Sequence();
			String seq_complementary = sequence.complementaryTrans(seq);
			// do seq RFLP
			boolean isIUPACenzyme = true;
			JudgeRFLP judgeRFLP = new JudgeRFLP(seq);
			judgeRFLP.setIUPACenzyme(isIUPACenzyme);
			if(judgeRFLP.isVarSeq()) {
				boolean isCut = judgeRFLP.isCanCut_dNTPs();
				resultList.add(isCut);
				List<List<String>> cutEnzymesList = judgeRFLP.getCutEnzymesList();
				cutEnzymesBatch.add(cutEnzymesList);
			}
			else {
				resultList.add(false);
			}
			// do seq_complementary RFLP
			JudgeRFLP judgeRFLP_complementary = new JudgeRFLP(seq_complementary);
			if(judgeRFLP_complementary.isVarSeq()) {
				boolean isCut = judgeRFLP_complementary.isCanCut_dNTPs();
				resultList_complementary.add(isCut);
				List<List<String>> cutEnzymesList_complementary = judgeRFLP_complementary.getCutEnzymesList();
				cutEnzymesBatch_complementary.add(cutEnzymesList_complementary);
			}
			else {
				resultList_complementary.add(false);
			}
		}
		/*for(int i=0;i<resultList.size();i++)
			System.out.println(resultList.get(i));*/
		//----------------------------------------------------------------------
		// Get elapsed time in milliseconds
		long elapsedTimeMillis = System.currentTimeMillis()-start;
		totalTime += elapsedTimeMillis;
		//System.out.println("Time: " + totalTime);
	}

	/**
	 * Get the result list for judge the variation if can be recognized by restriction enzymes.
	 * @return Result list for true or false.
	 */
	public List<Boolean> getIsCutResultList() {
		return resultList;
	}

	/**
	 * Get the complementary result list for judge the variation if can be recognized by restriction enzymes.
	 * @return Complementary result list for true or false.
	 */
	public List<Boolean> getIsCutResultListComplementary() {
		return resultList_complementary;
	}

	/**
	 * Get batch restriction enzymes.
	 * @return Batch restriction enzymes list.
	 */
	public List<List<List<String>>> getCutEnzymesBatch() {
		return cutEnzymesBatch;
	}

	/**
	 * Get batch complementary restriction enzymes.
	 * @return Batch complementary restriction enzymes list.
	 */
	public List<List<List<String>>> getCutEnzymesBatchComplementary() {
		return cutEnzymesBatch_complementary;
	}

	/**
	 * Get the total running time.
	 * @return The total running time.
	 */
	public long getTotalTime() {
		return totalTime;
	}

/*	public static void main(String args[]) throws Exception {
		List<String[]> seqList = new LinkedList<String[]>();
		String[] seq_data = new String[4];
		seq_data[0] = "8152226";
		seq_data[1] = "TATTCAAGTGCACGAGACCAATGAC";
		seq_data[2] = "[G/T]";
		seq_data[3] = "GGACCTCTGGTGAGGCCCTGGTGAG";
		seqList.add(seq_data);
		JudgeRFLPBatchThread judgeRFLPBatchThread = new JudgeRFLPBatchThread(seqList);
		judgeRFLPBatchThread.start();
	}*/
}