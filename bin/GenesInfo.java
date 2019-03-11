
public class GenesInfo implements Comparable<GenesInfo> {
	String geneName;
	double pValue;
	String geneSymbol;

	@Override
	public int compareTo(GenesInfo pSheetObject) {
		return pSheetObject.pValue < this.pValue ? 1 : (pSheetObject.pValue > this.pValue ? -1 : 0);
	}

	public String getGeneName() {
		return geneName;
	}

	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}

	public double getpValue() {
		return pValue;
	}

	public void setpValue(double pValue) {
		this.pValue = pValue;
	}

	public String getGeneSymbol() {
		return geneSymbol;
	}

	public void setGeneSymbol(String geneSymbol) {
		this.geneSymbol = geneSymbol;
	}
}
