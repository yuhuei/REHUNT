/*
 * Program name: RFLP.java
 * Date: 2016/09/30
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		RFLP analysis for a sequence.
 */

package bio.rehunt.rflp;

import java.io.*;
import java.util.*;

import bio.rehunt.seq.SeqProcess;

/**
 * RFLP analysis for a sequence.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 */
public class RFLP {
	// RFLPprocess
	private RFLPprocess process;
	private List<String> resEnzymeList1;	// 1.<ENZYME NAME>
	private List<String> resEnzymeList2;	// 2.<PROTOTYPE>
	private List<String> resEnzymeList3;	// 3.<MICROORGANISM>
	private List<String> resEnzymeList4;	// 4.<SOURCE>
	private List<String> resEnzymeList5;	// 5.<RECOGNITION SEQUENCE>
	private List<String> resEnzymeList6;	// 6.<METHYLATION SITE>
	private List<String> resEnzymeList7;	// 7.<COMMERCIAL AVAILABILITY>
	private List<String> resEnzymeList8;	// 8.<REFERENCES>
	private boolean isIUPACenzyme = false;	// judge whether get IUPAC enzymes or not
	private int enzymeLen_min = 4;
	private int enzymeLen_max = 10;
	private String seq;	// seq

	/**
	 * Constructor for initialize RFLP.
	 */
	public RFLP() {
		process = new RFLPprocess();
	}

	/**
	 * Set seq and filter unnecessary symbols.
	 * @param seq Sequence.
	 */
	// set the seq for search enzymes
	public void setSeq(String seq) {
		String upperSeq = seq.toUpperCase();
		StringBuffer strBuff_seq = new StringBuffer();
		for(int i=0;i<upperSeq.length();i++) {
			if(upperSeq.charAt(i)=='A' || upperSeq.charAt(i)=='T' || upperSeq.charAt(i)=='C' || upperSeq.charAt(i)=='G' ||
				upperSeq.charAt(i)=='M' || upperSeq.charAt(i)=='R' || upperSeq.charAt(i)=='W' || upperSeq.charAt(i)=='S' ||
				upperSeq.charAt(i)=='Y' || upperSeq.charAt(i)=='K' || upperSeq.charAt(i)=='V' || upperSeq.charAt(i)=='H' ||
				upperSeq.charAt(i)=='D' || upperSeq.charAt(i)=='B' || upperSeq.charAt(i)=='N') {
				strBuff_seq.append(upperSeq.charAt(i));
			}
		}
		this.seq = strBuff_seq.toString();
	}

	/**
	 * Get restriction enzymes from the sequence.
	 * @return If search successfully return true, else return false.
	 */
	public boolean getEnzymes() {
		// read restriction_enzyme file
		String strLine = "";
		InputStreamReader isr = null;
		BufferedReader bfr = null;
		try {
			isr = new InputStreamReader(getClass().getResourceAsStream("/REBASE/link_parsrefs.txt"));
			bfr = new BufferedReader(isr);
		}
		catch(Exception e) {
			System.out.println("Exception Message: " + e.getMessage());
			return false;
		}
		resEnzymeList1 = new LinkedList<String>();
		resEnzymeList2 = new LinkedList<String>();
		resEnzymeList3 = new LinkedList<String>();
		resEnzymeList4 = new LinkedList<String>();
		resEnzymeList5 = new LinkedList<String>();
		resEnzymeList6 = new LinkedList<String>();
		resEnzymeList7 = new LinkedList<String>();
		resEnzymeList8 = new LinkedList<String>();
		
		StringBuffer strBuff1 = new StringBuffer(15);
		StringBuffer strBuff2 = new StringBuffer(10);
		StringBuffer strBuff3 = new StringBuffer(55);
		StringBuffer strBuff4 = new StringBuffer(20);
		StringBuffer strBuff5 = new StringBuffer(35);
		StringBuffer strBuff6 = new StringBuffer(30);
		StringBuffer strBuff7 = new StringBuffer(20);
		StringBuffer strBuff8 = new StringBuffer(50);
		String enzymeSeq = null;

		int count = 0;
		try {
			while((strLine=bfr.readLine())!=null) {
				if(!strLine.equals("")) {
					if(strLine.charAt(0) != '<')
						continue;
					else {
						count++;
						// read data
						switch(strLine.charAt(1)) {
							case '1':
								strBuff1.append(strLine.substring(3));
								break;
							case '2':
								strBuff2.append(strLine.substring(3));
								break;
							case '3':
								strBuff3.append(strLine.substring(3));
								break;
							case '4':
								strBuff4.append(strLine.substring(3));
								break;
							case '5':
								strBuff5.append(strLine.substring(3));
								// process the enzyme recognition sequence
								enzymeSeq = process.removeResEnzymeSym(strBuff5.toString());
								break;
							case '6':
								strBuff6.append(strLine.substring(3));
								break;
							case '7':
								strBuff7.append(strLine.substring(3));
								break;
							case '8':
								strBuff8.append(strLine.substring(3));
								break;
							default:
								count--;
						}
					}
				}

				// process RFLP and clear strBuff
				if(count == 8) {
					// find the enzymes in the sequence
					if(!enzymeSeq.equals("")) {
						if(enzymeSeq.length() >= enzymeLen_min && enzymeSeq.length() <= enzymeLen_max) {
							if(process.haveRegionMatches(seq, enzymeSeq)) {
								resEnzymeList1.add(strBuff1.toString());
								resEnzymeList2.add(strBuff2.toString());
								resEnzymeList3.add(strBuff3.toString());
								resEnzymeList4.add(strBuff4.toString());
								resEnzymeList5.add(strBuff5.toString());
								resEnzymeList6.add(strBuff6.toString());
								resEnzymeList7.add(strBuff7.toString());
								resEnzymeList8.add(strBuff8.toString());
							}
							else if(isIUPACenzyme) {
								SeqProcess seqProcess = new SeqProcess();
								int allelePos = seqProcess.getFirstIUPACPos(enzymeSeq);
								if(allelePos != -1) {
									List<String> seqList = getAllIUPACseq(enzymeSeq);
									for(int i=0;i<seqList.size();i++) {
										//System.out.println("seqList: " + seqList.get(i).toString());
										if(process.haveRegionMatches(seq, seqList.get(i).toString())) {
											resEnzymeList1.add(strBuff1.toString());
											resEnzymeList2.add(strBuff2.toString());
											resEnzymeList3.add(strBuff3.toString());
											resEnzymeList4.add(strBuff4.toString());
											resEnzymeList5.add(strBuff5.toString());
											resEnzymeList6.add(strBuff6.toString());
											resEnzymeList7.add(strBuff7.toString());
											resEnzymeList8.add(strBuff8.toString());
											break;
										}
									}
								}
							}
						}
					}
					// clear strBuff
					strBuff1.delete(0, strBuff1.length());
					strBuff2.delete(0, strBuff2.length());
					strBuff3.delete(0, strBuff3.length());
					strBuff4.delete(0, strBuff4.length());
					strBuff5.delete(0, strBuff5.length());
					strBuff6.delete(0, strBuff6.length());
					strBuff7.delete(0, strBuff7.length());
					strBuff8.delete(0, strBuff8.length());
					count = 0;
				}
			}
			// close file
			bfr.close();
			isr.close();
		}
		catch(Exception e){
			return false;
		}
		return true;
	}

	/**
	 * Set the program if needs to find the IUPAC enzymes.
	 * @param isIUPACenzyme To set true is to find the IUPAC enzymes, and false is not.
	 */
	public void setIUPACenzyme(boolean isIUPACenzyme) {
		this.isIUPACenzyme = isIUPACenzyme;
	}

	/**
	 * Set minimum length of enzyme sequence for search.
	 * @param enzyme_len Restriction enzyme sequence length.
	 */
	public void setEnzymeLenMin(int enzyme_len) {
		enzymeLen_min = enzyme_len;
	}

	/**
	 * Set maximum length of enzyme sequence for search.
	 * @param enzyme_len Restriction enzyme sequence length.
	 */
	public void setEnzymeLenMax(int enzyme_len) {
		enzymeLen_max = enzyme_len;
	}

	/**
	 * Get all general sequences from IUPAC enzyme recognition sequence.
	 * @param IUPACseq The restriction enzyme sequence with IUPAC format.
	 * @return A sequence list from IUPAC enzyme recognition sequence.
	 */
	private List<String> getAllIUPACseq(String IUPACseq) {
		List<String> seqList = new LinkedList<String>();
		List<String> seqs = null;
		StringBuffer strBuff = new StringBuffer(IUPACseq);
		for(int i=0;i<IUPACseq.length();i++) {
			if(IUPACseq.charAt(i)=='M') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "A").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "C").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='R') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "A").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "G").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='W') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "A").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "T").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='S') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "C").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "G").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='Y') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "C").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "T").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='K') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "G").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "T").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='V') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "A").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "C").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "G").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='H') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "A").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "C").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "T").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='D') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "A").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "G").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "T").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='B') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "C").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "G").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "T").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
			else if(IUPACseq.charAt(i)=='N') {
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "A").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "C").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "G").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				seqs = getAllIUPACseq(strBuff.replace(i, i+1, "T").toString());
				for(int j=0;j<seqs.size();j++)
					seqList.add(seqs.get(j).toString());
				break;
			}
		}
		if(seqList.size() == 0)
		{
			seqList.add(IUPACseq);
			return seqList;
		}
		return seqList;
	}

	/**
	 * Get "ENZYME NAME" field list.
	 * @return Enzyme name list.
	 * @throws IOException IOException
	 */
	public List<String> getEnzymeNameList() throws IOException {
		return resEnzymeList1;
	}

	/**
	 * Get "PROTOTYPE" field list.
	 * @return Prototype list.
	 * @throws IOException IOException
	 */
	public List<String> getPrototypeList() throws IOException {
		return resEnzymeList2;
	}

	/**
	 * Get "MICROORGANISM" field list.
	 * @return Microorganism list.
	 * @throws IOException IOException
	 */
	public List<String> getMicroorganism() throws IOException {
		return resEnzymeList3;
	}

	/**
	 * Get "SOURCE" field list.
	 * @return Source list.
	 * @throws IOException IOException
	 */
	public List<String> getSource() throws IOException {
		return resEnzymeList4;
	}

	/**
	 * Get "RECOGNITION SEQUENCE" field list.
	 * @return Recognition seq list.
	 * @throws IOException IOException
	 */
	public List<String> getRecognition() throws IOException {
		return resEnzymeList5;
	}

	/**
	 * Get "METHYLATION SITE" field list.
	 * @return Methylation site list.
	 * @throws IOException IOException
	 */
	public List<String> getMethylation() throws IOException {
		return resEnzymeList6;
	}

	/**
	 * Get "COMMERCIAL AVAILABILITY" field list.
	 * @return Commercial availability list.
	 * @throws IOException IOException
	 */
	public List<String> getCommercial() throws IOException {
		return resEnzymeList7;
	}

	/**
	 * Get "REFERENCES" field list.
	 * @return References list.
	 * @throws IOException IOException
	 */
	public List<String> getReferences() throws IOException {
		return resEnzymeList8;
	}

/*	public static void main(String args[]) throws Exception
	{
		RFLP rflp = new RFLP();
		rflp.setSeq("RGA-TCY RGA-TC");
		rflp.setIUPACenzyme(true);
		rflp.setEnzymeLenMin(4);
		rflp.setEnzymeLenMax(10);
		if(rflp.getEnzymes())
		{
			List<String> enzymeNameList = rflp.getEnzymeNameList();
			System.out.println("enzymeNameList.size(): " + enzymeNameList.size());
//			for(int i=0;i<enzymeNameList.size();i++)
//				System.out.println(enzymeNameList.get(i));
//			List recognitionList = rflp.getRecognition();
//			for(int i=0;i<recognitionList.size();i++)
//				System.out.println(recognitionList.get(i));
		}
		else
			System.out.println("Error!");
	}*/
}