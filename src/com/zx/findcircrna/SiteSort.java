package com.zx.findcircrna;

public class SiteSort implements Comparable<SiteSort> {
	private Integer site;
	private String[] length;
	private int tieKey;

	public SiteSort(int site, String[] length) {
		this.site = site;
		this.length = length;
		this.tieKey = String.join("\t", length).hashCode();
	}

	public String[] getLength() {
		return length;
	}

	public Integer getSite() {
		return site;
	}

	@Override
	public int compareTo(SiteSort other) {
		int cmp = this.site - other.site;
		if (cmp != 0) return cmp;
		return Integer.compare(this.tieKey, other.tieKey);
	}
}
