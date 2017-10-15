/*
 * Program name: RFLPprocess.java
 * Date: 2016/09/30
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		RFLP process function.
 */

package bio.rehunt.rflp;

import java.util.*;

import bio.rehunt.algorithm.BM;

/**
 * RFLP process function.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 */
public class RFLPprocess {
	// Boyer-Moore
	private BM bm = new BM();

	/**
	 * Remove special symbols from a restriction enzyme_seq.
	 * @param enzymeSeq Restriction enzyme_seq sequence.
	 * @return The filtered restriction enzyme sequence.
	 */
	public String removeResEnzymeSym(String enzymeSeq) {
		StringBuffer resEnzyme = new StringBuffer();
		enzymeSeq = enzymeSeq.toUpperCase();
		for(int i=0;i<enzymeSeq.length();i++) {
			if(enzymeSeq.charAt(i) == 'A' || enzymeSeq.charAt(i) == 'C' || enzymeSeq.charAt(i) == 'G' ||
				enzymeSeq.charAt(i) == 'T' || enzymeSeq.charAt(i) == 'M' || enzymeSeq.charAt(i) == 'R' ||
				enzymeSeq.charAt(i) == 'W' || enzymeSeq.charAt(i) == 'S' || enzymeSeq.charAt(i) == 'Y' ||
				enzymeSeq.charAt(i) == 'K' || enzymeSeq.charAt(i) == 'V' || enzymeSeq.charAt(i) == 'H' ||
				enzymeSeq.charAt(i) == 'D' || enzymeSeq.charAt(i) == 'B' || enzymeSeq.charAt(i) == 'N')
				resEnzyme.append(enzymeSeq.charAt(i));
		}
		return resEnzyme.toString();
	}
	
	/**
	 * Count the sub sequence reappears in the sequence.
	 * @param seq The search sequence.
	 * @param subSeq The key sequence.
	 * @return The count of sub sequence reappears in the sequence.
	 */
	public int regionMatchCount(String seq, String subSeq) {
		int count = 0;
		bm.search(seq, subSeq);
		int pos = bm.search();
		if(pos != -1) {
			count++;
			while(bm.next() != -1)
				count++;
		}
		return count;
	}

	/**
	 * Judge if sequence contains the sub sequence.
	 * @param seq The search sequence.
	 * @param subSeq The key sequence.
	 * @return If the sequence contains the sub sequence then return true, else return flase.
	 */
	public boolean haveRegionMatches(String seq, String subSeq) {
		bm.search(seq, subSeq);
		int pos = bm.search();
		boolean find = false;
		if(pos != -1)
			find = true;
		return find;
	}

	/**
	 * Judge if sequence contains many sub sequence.
	 * @param seq The search sequence.
	 * @param subSeq The key sequence.
	 * @return If the sequence contains many sub sequence then return true, else return flase.
	 */
	public boolean haveSpecificityMatches(String seq, String subSeq) {
		bm.search(seq, subSeq);
		int pos = bm.search();
		boolean find = false;
		if(pos != -1) {
			if(bm.next() == -1)
				find = true;
		}
		return find;
	}

	/**
	 * Get the different elements between list1 and list2.
	 * @param list1 String list 1.
	 * @param list2 String list 2.
	 * @return Different elements between list1 and list2.
	 */
	public List<String> getDiffContent(List<String> list1, List<String> list2) {
		List<String> differentList = new LinkedList<String>();
		boolean[] list1_same = new boolean[list1.size()];
		boolean[] list2_same = new boolean[list2.size()];
		// initial
		for(int i=0;i<list1_same.length;i++)
			list1_same[i] = false;
		for(int j=0;j<list2_same.length;j++)
			list2_same[j] = false;
		// list_same = false represents the content does not appear to each other
		for(int i=0;i<list1.size();i++) {
			for(int j=0;j<list2.size();j++) {
				if(list2_same[j] == true)	// if list2_same[j] = true, not to compare
					continue;
	            if(list1.get(i).equals(list2.get(j))) {
	               list1_same[i] = true;
	               list2_same[j] = true;
	               break;
	            }
			}
		}
		for(int i=0;i<list1_same.length;i++) {
			if(list1_same[i] == false)
				differentList.add(list1.get(i));
		}
		for(int j=0;j<list2_same.length;j++) {
			if(list2_same[j] == false)
				differentList.add(list2.get(j));
		}
		return differentList;
	}

	/**
	 * Get the different elements between list1 and list2 and sort according to the alphabet order.
	 * @param list1 String list 1.
	 * @param list2 String list 2.
	 * @return Different elements that have sorted between list1 and list2.
	 */
	public List<String> getDiffOrder(List<String> list1, List<String> list2) {
		List<String> differentList = new LinkedList<String>();
		List<String> diffList1 = getDiff(list1, list2);
		List<String> diffList2 = getDiff(list2, list1);
		// sort by letters of the alphabet order
		Set<String> set = new TreeSet<String>();
		for(int i=0;i<diffList1.size();i++)
      		set.add(diffList1.get(i));
		for(int i=0;i<diffList2.size();i++)
      		set.add(diffList2.get(i));
		Iterator<String> iterator = set.iterator();
		while(iterator.hasNext()) {
			differentList.add(iterator.next().toString());
		}
		return differentList;
	}

	/**
	 * Get the different elements between main list and reference list and exclude reference list.
	 * @param mainList Main string list 1.
	 * @param refList Reference string list 2.
	 * @return Different elements between mainList and refList and exclude refList.
	 */
	public List<String> getDiff(List<String> mainList, List<String> refList) {
		boolean hit = false;
		List<String> diff = new LinkedList<String>();
		for(int i=0;i<mainList.size();i++) {
			for(int j=0;j<refList.size();j++) {
				if(mainList.get(i).equals(refList.get(j))) {
					hit = true;
					break;
				}
			}
			if(!hit)
				diff.add(mainList.get(i));
			hit = false;
		}
		return diff;
	}

	/**
	 * Get the same elements between main list and reference list.
	 * @param mainList Main string list 1.
	 * @param refList Reference string list 2.
	 * @return The same elements between main list and reference list.
	 */
	public List<String> getSame(List<String> mainList, List<String> refList) {
		boolean hit = false;
		List<String> same = new LinkedList<String>();
		for(int i=0;i<mainList.size();i++) {
			for(int j=0;j<refList.size();j++) {
				if(mainList.get(i).equals(refList.get(j))) {
					hit = true;
					break;
				}
			}
			if(hit)
				same.add(mainList.get(i));
			hit = false;
		}
		return same;
	}

	/**
	 * Segment data according to a symbol.
	 * @param data The data.
	 * @param symbol The symbol for segmenting data.
	 * @return The list for the segmented data.
	 */
	public List<String> segmentData(String data, String symbol) {
		List<String> dataList = new LinkedList<String>();
		int loc = 0;

		for(int i=0;i<data.length();i++) {
			for(int j=0;j<symbol.length();j++) {
				if(data.charAt(i) == symbol.charAt(j)) {
					if(loc != i)
						dataList.add(data.substring(loc, i).trim().toUpperCase());
					loc = i+1;
				}
			}
		}
		if(!data.substring(loc).equals(""))
			dataList.add(data.substring(loc).toUpperCase());
		return dataList;
	}

	/**
	 * Remove repeat elements.
	 * @param list String list.
	 * @return A string list without repeat elements.
	 */
	public List<String> removeRepeat(List<String> list) {
		List<String> noRepeat = new LinkedList<String>();
		if(list.size() == 0)
			return noRepeat;
		boolean isEqual = false;
		noRepeat.add(list.get(0));
		for(int i=1;i<list.size();i++) {
			for(int j=0;j<noRepeat.size();j++) {
				if(list.get(i).equals(noRepeat.get(j))) {
					isEqual = true;
					break;
				}
			}
			if(!isEqual)
				noRepeat.add(list.get(i));
			isEqual = false;
		}
		return noRepeat;
	}

/*	public static void main(String args[]) {
		// judge sequence if contains restriction-enzyme
		String sequence = "ACATGACTGAGGTCGTGAGACGCTG[C/G]CCCCACCATGAGCGTTGCTCTGATG";
		String reEnzyme = "TCGTGAGACGCT";
		RFLPprocess analysis = new RFLPprocess();
		System.out.println("haveRegionMatches: " + analysis.haveRegionMatches(sequence, reEnzyme));
		System.out.println("haveSpecificityMatches: " + analysis.haveSpecificityMatches(sequence, reEnzyme));
		// get the different enzyme name at list1 and list2
		List<String> list1 = new LinkedList<String>();
		List<String> list2 = new LinkedList<String>();
		list1.add("test");
		list1.add("seq");
		list1.add("reset");
		list1.add("about");
		list2.add("about");
		list2.add("union");
		list2.add("about");
		list2.add("seq");
		list2.add("visual");

		System.out.println("getDiffContent: " + analysis.getDiffContent(list1, list2));
		// process input data and return the data list
		String data = "rs111"+"\n\n\n,,,"+"rs222"+"\n"+"rs333"+"\r"+"rs888";
		char[] symbols = {',', '\n', '\r'};
		String symbol = new String(symbols);
		List<String> dataList = analysis.segmentData(data, symbol);
		for(int i=0;i<dataList.size();i++) {
			System.out.println(dataList.get(i));
		}
		// remove restriction enzyme special symbol
		String enzyme = "(10/15)GTacd^Cnnn";
		System.out.println(analysis.removeResEnzymeSym(enzyme));
		// get the different enzyme name that
		System.out.println("Diff: " + analysis.getDiff(list1, list2));
		System.out.println("DiffOrder: " + analysis.getDiffOrder(list1, list2));
		System.out.println("Same: " + analysis.getSame(list1, list2));
		// remove repeat
		System.out.println("NoRepeat: " + analysis.removeRepeat(list2));
	}*/
}