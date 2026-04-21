package com.zx.findcircrna;

public class IsStand {
	public String stand7(String infor) {
		return ((Integer.parseInt(infor) >> 6) & 1) == 0 ? "0" : "1";
	}
	public String stand5(String infor) {
		return ((Integer.parseInt(infor) >> 4) & 1) == 0 ? "0" : "1";
	}
}
