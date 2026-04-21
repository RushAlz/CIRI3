package com.zx.findcircrna;

public class SiteSort implements Comparable<SiteSort> {
	private Integer site;
	private String[] length;
	private int tieKey;

	public SiteSort(int site, String[] length) {
		this.site = site;
		this.length = length;
		this.tieKey = 0;
	}

	public SiteSort(int site, String[] length, int tieKey) {
		this.site = site;
		this.length = length;
		this.tieKey = tieKey;
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
		return other.tieKey - this.tieKey;
	}
}
