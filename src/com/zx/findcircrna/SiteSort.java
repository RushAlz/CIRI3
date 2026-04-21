package com.zx.findcircrna;

public class SiteSort implements Comparable<SiteSort> {
	private Integer site;
	private String[] length;

	public SiteSort(int site, String[] length) {
		this.site = site;
		this.length = length;
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
		cmp = Integer.parseInt(this.length[0]) - Integer.parseInt(other.length[0]);
		if (cmp != 0) return cmp;
		return Integer.parseInt(this.length[1]) - Integer.parseInt(other.length[1]);
	}
}
