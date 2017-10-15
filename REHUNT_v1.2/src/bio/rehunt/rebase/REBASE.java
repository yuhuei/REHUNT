/*
 * Program name: REBASE.java
 * Date: 2016/09/30
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 *		Search the restriction enzyme related information using enzyme name from REBASE link_parsrefs.txt.
 */

package bio.rehunt.rebase;

import java.io.*;
import java.util.*;

import bio.rehunt.rflp.RFLPprocess;

/**
 * Search the restriction enzyme related information using enzyme name from REBASE link_parsrefs.txt.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 */
public class REBASE {
	// reference
	private String authors = null;
	private String year = null;
	private String journal = null;
	private String volume = null;
	private String page = null;
	private String unpublish = null;

	// file read
	private InputStreamReader isr = null;
	private BufferedReader bfr = null;

	/**
	 * Get file reader for REBASE link_parsrefs.txt.
	 * @return If success return true, else return false.
	 */
	public boolean getRebase() {
		try {
			isr = new InputStreamReader(getClass().getResourceAsStream("/REBASE/link_parsrefs.txt"));
			bfr = new BufferedReader(isr);
		}
		catch(Exception e) {
			return false;
		}
		return true;
	}

	/**
	 * Search enzyme information using enzyme name.
	 * @param enzymeName Restriction enzyme name.
	 * @return Enzyme information.
	 * @throws IOException IOException
	 */
	public String[] searchEnzyme(String enzymeName) throws IOException {
		boolean isSuccessRebase = getRebase();
		if(!isSuccessRebase)
			return null;
		String[] data = new String[8];
		String strLine = "";
		int count = 0;
		while((strLine=bfr.readLine()) != null) {
			if(!strLine.equals("")) {
				if(strLine.charAt(0) != '<')
					continue;
				else {
					count++;
					// read data
					switch(strLine.charAt(1)) {
						case '1':
							data[0] = strLine.substring(3);
							break;
						case '2':
							data[1] = strLine.substring(3);
							break;
						case '3':
							data[2] = strLine.substring(3);
							break;
						case '4':
							data[3] = strLine.substring(3);
							break;
						case '5':
							data[4] = strLine.substring(3);
							break;
						case '6':
							data[5] = strLine.substring(3);
							break;
						case '7':
							data[6] = strLine.substring(3);
							break;
						case '8':
							data[7] = strLine.substring(3);
							break;
						default:
							count--;
					}
					if(enzymeName.equals(data[0]) && count==8)
						break;
				}
			}
			else {	// clear
				data[0] = "";
				data[1] = "";
				data[2] = "";
				data[3] = "";
				data[4] = "";
				data[5] = "";
				data[6] = "";
				data[7] = "";
				count = 0;
			}
		}
		setClose();
		if(data[0] == null)
			return null;
		return data;
	}

	/**
	 * Search enzyme information using enzyme recognition sequence.
	 * @param recognition Restriction enzyme recognition sequence.
	 * @return Enzyme information list.
	 * @throws IOException IOException
	 */
	public List<String[]> searchRecognition(String recognition) throws IOException {
		boolean isSuccessRebase = getRebase();
		if(!isSuccessRebase)
			return null;
		RFLPprocess process = new RFLPprocess();
		List<String[]> dataList = new LinkedList<String[]>();
		String[] data = new String[8];
		String strLine = "";
		int count = 0;
		while((strLine=bfr.readLine()) != null) {
			if(!strLine.equals("")) {
				if(strLine.charAt(0) != '<')
					continue;
				else {
					count++;
					// read data
					switch(strLine.charAt(1)) {
						case '1':
							data[0] = strLine.substring(3);
							break;
						case '2':
							data[1] = strLine.substring(3);
							break;
						case '3':
							data[2] = strLine.substring(3);
							break;
						case '4':
							data[3] = strLine.substring(3);
							break;
						case '5':
							data[4] = strLine.substring(3);
							break;
						case '6':
							data[5] = strLine.substring(3);
							break;
						case '7':
							data[6] = strLine.substring(3);
							break;
						case '8':
							data[7] = strLine.substring(3);
							break;
						default:
							count--;
					}
					if(data[4] != null && count==8) {
						if(recognition.equals(process.removeResEnzymeSym(data[4]))) {
							String[] dataCopy = data.clone();
							dataList.add(dataCopy);
						}
					}
				}
			}
			else {	// clear
				data[0] = "";
				data[1] = "";
				data[2] = "";
				data[3] = "";
				data[4] = "";
				data[5] = "";
				data[6] = "";
				data[7] = "";
				count = 0;
			}
		}
		setClose();
		if(dataList.size() == 0)
			return null;
		return dataList;
	}

	/**
	 * Get "REFERENCES" string.
	 * @param index Index of restriction enzyme reference.
	 * @return Enzyme reference string.
	 * @throws IOException IOException
	 */
	public String getReferenceStr(String index) throws IOException {
		boolean isSuccessRebase = getRebase();
		if(!isSuccessRebase)
			return null;
		String data = "";
		String strLine = "";
		boolean start = false;
		while((strLine=bfr.readLine())!=null) {
			if(!start) {
				if(!strLine.equals("References:"))
					continue;
				else {
					start = true;
					continue;
				}
			}
			else {
				if(!strLine.equals("")) {
					//System.out.println(strLine);
					if(strLine.substring(0, strLine.indexOf('.')).equals(index)) {
						int authorsPos = strLine.indexOf("<AUTHORS>");
						int jPos = strLine.indexOf("<J>");
						int yearPos = 0;
						int journalPos = 0;
						int volumePos = 0;
						int pagePos = 0;
						int uPos = 0;
						if(jPos != -1) {
							yearPos = strLine.indexOf("<YEAR>");
							journalPos = strLine.indexOf("<JOURNAL>");
							volumePos = strLine.indexOf("<VOLUME>");
							pagePos = strLine.indexOf("<PAGES>");
							authors = strLine.substring(authorsPos+9, jPos).trim();
							if(yearPos != -1) {
								if(strLine.indexOf("<", yearPos+6) != -1)
									year = strLine.substring(yearPos+6, strLine.indexOf("<", yearPos+6)).trim();
								else
									year = strLine.substring(yearPos+6).trim();
							}
							if(journalPos != -1) {
								if(strLine.indexOf("<", journalPos+9) != -1)
									journal = strLine.substring(journalPos+9, strLine.indexOf("<", journalPos+9)).trim();
								else
									journal = strLine.substring(journalPos+9).trim();
							}
							if(volumePos != -1) {
								if(strLine.indexOf("<", volumePos+8) != -1)
									volume = strLine.substring(volumePos+8, strLine.indexOf("<", volumePos+8)).trim();
								else
									volume = strLine.substring(volumePos+8).trim();
							}
							if(volumePos != -1)
								page = strLine.substring(pagePos+7).trim();
							data = authors + ", (" +  year + ") " + journal + ", vol. " + volume + ", pp. " + page + ".";
						}
						else {
							uPos = strLine.indexOf("<U>");
							authors = strLine.substring(authorsPos+9, uPos).trim();
							unpublish = strLine.substring(uPos+3).trim();
							data = authors + ", " +  unpublish;
						}
					}
				}
			}
		}
		setClose();
		return data;
	}

	/**
	 * Get "ENZYME NAME" string.
	 * @param data Enzyme information array.
	 * @return Enzyme name string.
	 */
	public String getEnzymeName(String[] data) {
		return data[0];
	}

	/**
	 * Get "PROTOTYPE" string.
	 * @param data Enzyme information array.
	 * @return Prototype string.
	 */
	public String getPrototype(String[] data) {
		return data[1];
	}

	/**
	 * Get "MICROORGANISM" string.
	 * @param data Enzyme information array.
	 * @return Microorganism string.
	 */
	public String getMicroorganism(String[] data) {
		return data[2];
	}

	/**
	 * Get "SOURCE" string.
	 * @param data Enzyme information array.
	 * @return Source string.
	 */
	public String getSource(String[] data) {
		return data[3];
	}

	/**
	 * Get "RECOGNITION SEQUENCE" string.
	 * @param data Enzyme information array.
	 * @return Recognition sequence string.
	 */
	public String getRecognition(String[] data) {
		return data[4];
	}

	/**
	 * Get "METHYLATION SITE" string.
	 * @param data Enzyme information array.
	 * @return Methylation site string.
	 */
	public String getMethylation(String[] data) {
		return data[5];
	}

	/**
	 * Get "COMMERCIAL AVAILABILITY" string.
	 * @param data Enzyme information array.
	 * @return Commercial availability string.
	 */
	public String getCommercial(String[] data) {
		return data[6];
	}

	/**
	 * Get "REFERENCES" indices.
	 * @param data Enzyme information array.
	 * @return References indices string with distinguished symbol ','.
	 */
	public String getReferences(String[] data) {
		return data[7];
	}

	/**
	 * Close REBASE database file.
	 * @throws IOException IOException
	 */
	// close file
	public void setClose() throws IOException {
		// close file
		if(bfr != null)
			bfr.close();
		if(isr != null)
			isr.close();
	}

/*	public static void main(String args[]) throws IOException {
//		// search from enzyme name
//		REBASE rebase = new REBASE();
//		String data[] = rebase.searchEnzyme("BsiHKAI");
//		if(data!=null) {
//			for(int i=0;i<data.length;i++)
//				System.out.println(data[i]);
//		}
//		else
//			System.out.println("Not to found...");
//
//		String enzymeName = rebase.getEnzymeName(data);
//		String prototype = rebase.getPrototype(data);
//		String microorganism = rebase.getMicroorganism(data);
//		String source = rebase.getSource(data);
//		String recognition = rebase.getRecognition(data);
//		String methylation = rebase.getMethylation(data);
//		String commercial = rebase.getCommercial(data);
//		String references = rebase.getReferences(data);
//		System.out.println(enzymeName);
//		System.out.println(prototype);
//		System.out.println(microorganism);
//		System.out.println(source);
//		System.out.println(recognition);
//		System.out.println("methylation:" + methylation);
//		System.out.println("commercial: " + commercial);
//		System.out.println(references);
//		if(references != null) {
//			String[] referencesIndex = references.split(",");
//			for(int i=0;i<referencesIndex.length;i++)
//				System.out.println("[" + (i+1) + "] " + rebase.getReferenceStr(referencesIndex[i]));
//		}
//
//		data = rebase.searchEnzyme("AaaI");
//		if(data!=null) {
//			for(int i=0;i<data.length;i++)
//				System.out.println(data[i]);
//		}
//		else
//			System.out.println("Not to found...");

		// search from enzyme recognition
		REBASE rebase = new REBASE();
		List<String[]> dataList = rebase.searchRecognition("GCWGC");
		if(dataList != null) {
			for(int i=0;i<dataList.size();i++) {
				String[] data = dataList.get(i);
				String enzymeName = rebase.getEnzymeName(data);
				String prototype = rebase.getPrototype(data);
				String microorganism = rebase.getMicroorganism(data);
				String source = rebase.getSource(data);
				String recognition = rebase.getRecognition(data);
				String methylation = rebase.getMethylation(data);
				String commercial = rebase.getCommercial(data);
				String references = rebase.getReferences(data);
				System.out.println(enzymeName);
				System.out.println(prototype);
				System.out.println(microorganism);
				System.out.println(source);
				System.out.println(recognition);
				System.out.println("methylation:" + methylation);
				System.out.println("commercial: " + commercial);
				System.out.println(references);
				if(references != null) {
					String[] referencesIndex = references.split(",");
					for(int j=0;j<referencesIndex.length;j++)
						System.out.println("[" + (j+1) + "] " + rebase.getReferenceStr(referencesIndex[j]));
				}
			}
		}
		else
			System.out.println("Not to found...");

		dataList = rebase.searchRecognition("GCWGC");
		if(dataList != null) {
			for(int i=0;i<dataList.size();i++) {
				String[] data = dataList.get(i);
				for(int j=0;j<data.length;j++)
					System.out.println(data[j]);
			}
		}
		else
			System.out.println("Not to found...");
	}*/
}