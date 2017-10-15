/*
 * Program name: BM.java
 * Date: 2016/09/30
 * Author: Yu-Huei Cheng
 * E-mail: yuhuei.cheng@gmail.com
 *
 * Function:
 * 		Boyer-Moore algorithm for string search.
 */

package bio.rehunt.algorithm;

import java.io.*;

/**
 * Boyer-Moore algorithm for string search.
 * @author Yu-Huei Cheng
 * @version REHUNT v1.2
 * @since JDK1.8.0
 */
public class BM {
 	private String text;	// total text
 	private String key;	// key string
 	private int[] skip;	// skip table
 	private int pos;	// current search position
 	
	/**
	 * Constructor for initialize boyer-moore algorithm.
	 */
 	public BM() {
 		text = "";
 		key = "";
 		skip = new int[256];
 		pos = 0;
 	}
 	
	/**
	 * Constructor for initialize boyer-moore algorithm.
	 * @param text The search text.
	 * @param key The search key.
	 */
 	public BM(String text, String key) {
 		this.text = text;
 		this.key = key;
 		skip = new int[256];
 		pos = 0;
 		setSkip();	// set skip table
 	}
 	
	/**
	 * Search key from text. The text and key must be set firstly.
	 * @return The found key position on text.
	 */
 	public int search() {
 		pos = search(0, text, key);
 		return getFindPos();
 	}
 	
	/**
	 * Search key from text. The key must be set firstly.
	 * @param text The search text.
	 * @return The found key position on text.
	 */
 	public int search(String text) {
 		this.text = text;
 		pos = search(0, text, key);
 		return getFindPos();
 	}
 	
	/**
	 * Search key from text.
	 * @param text The search text.
	 * @param key The search key.
	 * @return The found key position on text.
	 */
 	public int search(String text, String key) {
 		this.text = text;
 		this.key = key;
 		setSkip();
 		pos = search(0, text, key);
 		return getFindPos();
 	}
 	
	/**
	 * Search key from text.
	 * @param pos The start position for search.
	 * @param text The search text.
	 * @param key The search key.
	 * @return The found position of key on text.
	 */
 	public int search(int pos, String text, String key) {
 		if(key.equals(""))
 			return -1;
		while(pos+key.length() <= text.length()) {
			String tmp = text.substring(pos, pos+key.length());
			if(tmp.equals(key)) {
				pos += key.length();
				return pos;
			}
			pos += skip[text.charAt(pos+key.length()-1)];
		}
		return -1;
 	}
 	
	/**
	 * Set skip table.
	 */
 	private void setSkip() {
 		for(int i=0;i<256;i++)	// initial skip length is the key length
 			skip[i] = key.length();
 		for(int i=0;i<key.length()-1;i++)	// adjust the skip
 			skip[key.charAt(i)] = key.length()-i-1;
 	}
 	
	/**
	 * Search next key.
	 * @return The found next position of key on text.
	 */
	public int next() {
		pos = search(pos, text, key);
		return getFindPos();
	}
	
	/**
	 * Set search text.
	 * @param text The search text.
	 */
	public void setText(String text) {
		this.text =	text;
	}
	
	/**
	 * Set search key.
	 * @param key The search key.
	 */
	public void setKey(String key) {
		this.key = key;
		setSkip();
	}
	
	/**
	 * Get found position of key on text.
	 * @return The found position of key on text.
	 */
	private int getFindPos() {
		int find_pos;	// find position
		if(pos != -1)
			find_pos = pos - key.length();
		else
			find_pos = -1;
		return find_pos;
	}
	
	/**
	 * Reset the start position for search.
	 */
	public void reset() {
		pos = 0;
	}
	
	public static void main(String args[]) throws Exception {
		// read test_file
		String test_file = "test" + File.separator + "test_file.fasta";
		FileReader fr = new FileReader(test_file);
		BufferedReader bfr = new BufferedReader(fr);
		String str = "";
		StringBuffer strBuff = new StringBuffer();
		while((str=bfr.readLine()) != null) {
			strBuff.append(str);
			strBuff.append("\n");
		}
		bfr.close();
		String text = strBuff.toString();
		String key = "CCCCC";
		System.out.println("len: " + text.length());
		//----------------------------------------------------------------------
		// Examples for brute force
		//----------------------------------------------------------------------
		int pos_bf = 0;
		for(int i=0;i<text.length();i++) {
			if(pos_bf+key.length() < text.length()) {
				String tmp = text.substring(pos_bf, pos_bf+key.length());
				if(tmp.equals(key)) {
					System.out.println("The key is found on: " + pos_bf);
				}
				pos_bf++;
			}
		}
		
		//----------------------------------------------------------------------
		// Examples for boyer moore
		//----------------------------------------------------------------------
		// example 1
		BM bm = new BM(text, key);
		int pos_bm = 0;
		while((pos_bm = bm.next()) != -1) {
			System.out.println("The key is found on: " + pos_bm);
		}
		
		/*// example 2
		BM bm = new BM();
		bm.setText(text);
		bm.setKey(key);
		int pos_bm = 0;
		while((pos_bm = bm.next()) != -1) {
			System.out.println("The key is found on: " + pos_bm);
		}*/
	}
}