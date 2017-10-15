/*
 * Program name: EnzymeClassification.java
 * Date: 2016/09/30
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		Enzymes classfication.
 */

package bio.rehunt.rebase;

import java.io.*;
import java.util.*;

import bio.rehunt.seq.SeqProcess;
import bio.rehunt.rflp.RFLPprocess;

/**
 * Enzymes classfication.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 */
public class EnzymeClassification {
	private REBASE rebase = null;
	private SeqProcess seqProcess = null;
	private RFLPprocess rflpProcess = null;
	private List<String> enzymeSeq_list = null;
	private List<String> enzymeNames_commercial = null;
	private List<String> enzymeNames_nonCommercial = null;
	// commercial and nonCommercial for IUPAC and no IUPAC
	private List<String> enzymeSeq_iupac = null;
	private List<String> enzymeSeq_noiupac = null;
	private List<String> enzymeSeq_iupac_commercial = null;
	private List<String> enzymeSeq_noiupac_commercial = null;
	private List<String> enzymeSeq_iupac_nonCommercial = null;
	private List<String> enzymeSeq_noiupac_nonCommercial = null;
	private List<String> enzymeNames_iupac_commercial = null;
	private List<String> enzymeNames_noiupac_commercial = null;
	private List<String> enzymeNames_iupac_nonCommercial = null;
	private List<String> enzymeNames_noiupac_nonCommercial = null;

	/**
	 * Constructor for initialize EnzymeClassification.
	 */
	public EnzymeClassification() {
		this.rebase = new REBASE();
		this.seqProcess = new SeqProcess();
		this.rflpProcess = new RFLPprocess();
	}

	/**
	 * Classfy restriction enzymes with commercial, non-commercial and IUPAC format recognition sequence.
	 * @param enzymeName_list Restriction enzyme name list.
	 * @return If classfy successfully then return true, else return false.
	 */
	public boolean classifyCommercialEnzymes(List<String> enzymeName_list) {
		enzymeSeq_list = new LinkedList<String>();
		enzymeNames_commercial = new LinkedList<String>();
		enzymeNames_nonCommercial = new LinkedList<String>();
		// commercial and nonCommercial for IUPAC and no IUPAC
		enzymeSeq_iupac = new LinkedList<String>();
		enzymeSeq_noiupac = new LinkedList<String>();
		enzymeSeq_iupac_commercial = new LinkedList<String>();
		enzymeSeq_noiupac_commercial = new LinkedList<String>();
		enzymeSeq_iupac_nonCommercial = new LinkedList<String>();
		enzymeSeq_noiupac_nonCommercial = new LinkedList<String>();
		enzymeNames_iupac_commercial = new LinkedList<String>();
		enzymeNames_noiupac_commercial = new LinkedList<String>();
		enzymeNames_iupac_nonCommercial = new LinkedList<String>();
		enzymeNames_noiupac_nonCommercial = new LinkedList<String>();
		try {
			for(int i=0;i<enzymeName_list.size();i++) {
				String enzymeName = enzymeName_list.get(i);
				String data[] = rebase.searchEnzyme(enzymeName);
				//String enzymeName = rebase.getEnzymeName(data);
				String recognition = rebase.getRecognition(data);
				String commercial = rebase.getCommercial(data);
				String enzymeSeq = rflpProcess.removeResEnzymeSym(recognition);
				enzymeSeq_list.add(enzymeSeq);
				//System.out.println("COMMERCIAL AVAILABILITY: " + commercial);
				if(!commercial.equals("")) {
					//System.out.println("Commercial");
					enzymeNames_commercial.add(enzymeName);
					if(seqProcess.getFirstIUPACPos(enzymeSeq) != -1) {
						enzymeSeq_iupac.add(enzymeSeq);
						enzymeSeq_iupac_commercial.add(enzymeSeq);
						enzymeNames_iupac_commercial.add(enzymeName);
					}
					else {
						enzymeSeq_noiupac.add(enzymeSeq);
						enzymeSeq_noiupac_commercial.add(enzymeSeq);
						enzymeNames_noiupac_commercial.add(enzymeName);
					}
				}
				else {
					//System.out.println("Non-Commercial");
					enzymeNames_nonCommercial.add(enzymeName);
					if(seqProcess.getFirstIUPACPos(enzymeSeq) != -1) {
						enzymeSeq_iupac.add(enzymeSeq);
						enzymeSeq_iupac_nonCommercial.add(enzymeSeq);
						enzymeNames_iupac_nonCommercial.add(enzymeName);
					}
					else {
						enzymeSeq_noiupac.add(enzymeSeq);
						enzymeSeq_noiupac_nonCommercial.add(enzymeSeq);
						enzymeNames_noiupac_nonCommercial.add(enzymeName);
					}
				}
			}
		}
		catch(IOException e) {
			System.out.println(e);
			return false;
		}
		return true;
	}

	/**
	 * Get commercial restriction enzyme name list.
	 * @return Commercial restriction enzyme name list.
	 */
	public List<String> getEnzymeNamesCommercial() {
		return enzymeNames_commercial;
	}

	/**
	 * Get non-commercial restriction enzyme name list.
	 * @return Non-commercial restriction enzyme name list.
	 */
	public List<String> getEnzymeNamesNonCommercial() {
		return enzymeNames_nonCommercial;
	}

	/**
	 * Get restriction enzyme recognition sequence list without repeat.
	 * @return Restriction enzyme recognition sequence list without repeat.
	 */
	public List<String> getEnzymeSeqList() {
		enzymeSeq_list = rflpProcess.removeRepeat(enzymeSeq_list);
		return enzymeSeq_list;
	}

	/**
	 * Get IUPAC restriction enzyme recognition sequence list without repeat.
	 * @return IUPAC restriction enzyme recognition sequence list without repeat.
	 */
	public List<String> getEnzymeSeqIupac() {
		enzymeSeq_iupac = rflpProcess.removeRepeat(enzymeSeq_iupac);
		return enzymeSeq_iupac;
	}

	/**
	 * Get non-IUPAC restriction enzyme recognition sequence list without repeat.
	 * @return Non-IUPAC restriction enzyme recognition sequence list without repeat.
	 */
	public List<String> getEnzymeSeqNoIupac() {
		enzymeSeq_noiupac = rflpProcess.removeRepeat(enzymeSeq_noiupac);
		return enzymeSeq_noiupac;
	}

	/**
	 * Get commercial IUPAC restriction enzyme recognition sequence list without repeat.
	 * @return Commercial IUPAC restriction enzyme recognition sequence list without repeat.
	 */
	public List<String> getEnzymeSeqIupacCommercial() {
		enzymeSeq_iupac_commercial = rflpProcess.removeRepeat(enzymeSeq_iupac_commercial);
		return enzymeSeq_iupac_commercial;
	}

	/**
	 * Get commercial non-IUPAC restriction enzyme recognition sequence list without repeat.
	 * @return Commercial non-IUPAC restriction enzyme recognition sequence list without repeat.
	 */
	public List<String> getEnzymeSeqNoIupacCommercial() {
		enzymeSeq_noiupac_commercial = rflpProcess.removeRepeat(enzymeSeq_noiupac_commercial);
		return enzymeSeq_noiupac_commercial;
	}

	/**
	 * Get Non-commercial IUPAC restriction enzyme recognition sequence list without repeat.
	 * @return Non-commercial IUPAC restriction enzyme recognition sequence list without repeat.
	 */
	public List<String> getEnzymeSeqIupacNonCommercial() {
		enzymeSeq_iupac_nonCommercial = rflpProcess.removeRepeat(enzymeSeq_iupac_nonCommercial);
		return enzymeSeq_iupac_nonCommercial;
	}

	/**
	 * Get non-commercial non-IUPAC restriction enzyme recognition sequence list without repeat.
	 * @return Non-commercial non-IUPAC restriction enzyme recognition sequence list without repeat.
	 */
	public List<String> getEnzymeSeqNoIupacNonCommercial() {
		enzymeSeq_noiupac_nonCommercial = rflpProcess.removeRepeat(enzymeSeq_noiupac_nonCommercial);
		return enzymeSeq_noiupac_nonCommercial;
	}

	/**
	 * Get commercial IUPAC restriction enzyme name list without repeat.
	 * @return Commercial IUPAC restriction enzyme name list without repeat.
	 */
	public List<String> getEnzymeNamesIupacCommercial() {
		return enzymeNames_iupac_commercial;
	}

	/**
	 * Get commercial non-IUPAC restriction enzyme name list without repeat.
	 * @return Commercial non-IUPAC restriction enzyme name list without repeat.
	 */
	public List<String> getEnzymeNamesNoIupacCommercial() {
		return enzymeNames_noiupac_commercial;
	}

	/**
	 * Get non-commercial IUPAC restriction enzyme name list without repeat.
	 * @return Non-commercial IUPAC restriction enzyme name list without repeat.
	 */
	public List<String> getEnzymeNamesIupacNonCommercial() {
		return enzymeNames_iupac_nonCommercial;
	}

	/**
	 * Get non-commercial non-IUPAC restriction enzyme name list without repeat.
	 * @return Non-commercial non-IUPAC restriction enzyme name list without repeat.
	 */
	public List<String> getEnzymeNamesNoIupacNonCommercial() {
		return enzymeNames_noiupac_nonCommercial;
	}

/*	public static void main(String args[]) throws IOException {
		List<String> enzymeName_list = new LinkedList<String>();
		enzymeName_list.add("AagI");
		enzymeName_list.add("CviJI");
		enzymeName_list.add("CviKI-1");

		EnzymeClassification enzymeClassification = new EnzymeClassification();
		if(enzymeClassification.classifyCommercialEnzymes(enzymeName_list)) {
			// get enzymes for commercial and nonCommercial
			List<String> enzymeSeq_list = enzymeClassification.getEnzymeSeqList();
			List<String> enzymeNames_commercial = enzymeClassification.getEnzymeNamesCommercial();
			List<String> enzymeNames_nonCommercial = enzymeClassification.getEnzymeNamesNonCommercial();
			for(int i=0;i<enzymeSeq_list.size();i++)
				System.out.println("Enzyme seq: " + enzymeSeq_list.get(i));
			for(int i=0;i<enzymeNames_commercial.size();i++)
				System.out.println("Enzyme commercial: " + enzymeNames_commercial.get(i));
			for(int i=0;i<enzymeNames_nonCommercial.size();i++)
				System.out.println("Enzyme noncommercial: " + enzymeNames_nonCommercial.get(i));
			// get enzymes for commercial and nonCommercial of IUPAC or non-IUPAC
			List<String> enzymeSeq_iupac = enzymeClassification.getEnzymeSeqIupac();
			List<String> enzymeSeq_noiupac = enzymeClassification.getEnzymeSeqNoIupac();
			List<String> enzymeSeq_iupac_commercial = enzymeClassification.getEnzymeSeqIupacCommercial();
			List<String> enzymeSeq_noiupac_commercial = enzymeClassification.getEnzymeSeqNoIupacCommercial();
			List<String> enzymeSeq_iupac_nonCommercial = enzymeClassification.getEnzymeSeqIupacNonCommercial();
			List<String> enzymeSeq_noiupac_nonCommercial = enzymeClassification.getEnzymeSeqNoIupacNonCommercial();
			List<String> enzymeNames_iupac_commercial = enzymeClassification.getEnzymeNamesIupacCommercial();
			List<String> enzymeNames_noiupac_commercial = enzymeClassification.getEnzymeNamesNoIupacCommercial();
			List<String> enzymeNames_iupac_nonCommercial = enzymeClassification.getEnzymeNamesIupacNonCommercial();
			List<String> enzymeNames_noiupac_nonCommercial = enzymeClassification.getEnzymeNamesNoIupacNonCommercial();
			for(int i=0;i<enzymeSeq_iupac.size();i++)
				System.out.println("IUPAC seq: " + enzymeSeq_iupac.get(i));
			for(int i=0;i<enzymeSeq_noiupac.size();i++)
				System.out.println("NO IUPAC seq: " + enzymeSeq_noiupac.get(i));
			for(int i=0;i<enzymeSeq_iupac_commercial.size();i++)
				System.out.println("IUPAC commercial seq: " + enzymeSeq_iupac_commercial.get(i));
			for(int i=0;i<enzymeSeq_noiupac_commercial.size();i++)
				System.out.println("NO IUPAC commercial seq: " + enzymeSeq_noiupac_commercial.get(i));
			for(int i=0;i<enzymeSeq_iupac_nonCommercial.size();i++)
				System.out.println("IUPAC nonCommercial seq: " + enzymeSeq_iupac_nonCommercial.get(i));
			for(int i=0;i<enzymeSeq_noiupac_nonCommercial.size();i++)
				System.out.println("No IUPAC nonCommercial seq: " + enzymeSeq_noiupac_nonCommercial.get(i));
			for(int i=0;i<enzymeNames_iupac_commercial.size();i++)
				System.out.println("IUPAC commercial: " + enzymeNames_iupac_commercial.get(i));
			for(int i=0;i<enzymeNames_noiupac_commercial.size();i++)
				System.out.println("NO IUPAC commercial: " + enzymeNames_noiupac_commercial.get(i));
			for(int i=0;i<enzymeNames_iupac_nonCommercial.size();i++)
				System.out.println("IUPAC nonCommercial: " + enzymeNames_iupac_nonCommercial.get(i));
			for(int i=0;i<enzymeNames_noiupac_nonCommercial.size();i++)
				System.out.println("No IUPAC nonCommercial: " + enzymeNames_noiupac_nonCommercial.get(i));
		}
	}*/
}