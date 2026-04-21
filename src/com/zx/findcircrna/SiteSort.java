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
		// Primary: compare by site.
		int c = this.site - other.site;
		if (c != 0) return c;
		// Tiebreak: compare the coordinate array lexicographically.
		// Without this, two SiteSort entries that share the same site end up
		// in insertion order, which depends on the HashSet iteration order
		// that produced them — non-deterministic across JVM runs. Scan2's
		// BSJ-detection code (IsBSJScan2Star / IsBSJScan2Intron / etc.)
		// walks the sorted list and returns the FIRST matching circRNA, so
		// ties resolving differently change which circRNA a read is
		// attributed to. Making the comparator a total order eliminates
		// that source of run-to-run drift.
		String[] a = this.length, b = other.length;
		if (a == null && b == null) return 0;
		if (a == null) return -1;
		if (b == null) return  1;
		int n = Math.min(a.length, b.length);
		for (int i = 0; i < n; i++) {
			c = a[i].compareTo(b[i]);
			if (c != 0) return c;
		}
		return a.length - b.length;
	}
}
