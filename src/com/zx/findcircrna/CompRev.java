package com.zx.findcircrna;

public class CompRev {

	private static final char[] COMP = new char[128];
	static {
		COMP['A'] = 'T'; COMP['T'] = 'A';
		COMP['G'] = 'C'; COMP['C'] = 'G';
		COMP['a'] = 't'; COMP['t'] = 'a';
		COMP['g'] = 'c'; COMP['c'] = 'g';
		COMP['N'] = 'N'; COMP['n'] = 'n';
	}

	public String compRev(String seq) {
		int n = seq.length();
		char[] out = new char[n];
		for (int i = 0; i < n; i++) {
			char c = seq.charAt(n - 1 - i);
			out[i] = (c < 128 && COMP[c] != 0) ? COMP[c] : c;
		}
		return new String(out);
	}
}
