import java.util.List;

public class StringListSorting implements Comparable<StringListSorting> {
	List<String> geneName;
	double pValue;;

	@Override
	public int compareTo(StringListSorting StringListSortingbject) {
		return StringListSortingbject.pValue < this.pValue ? 1 : (StringListSortingbject.pValue > this.pValue ? -1 : 0);
	}

	public List<String> getGeneName() {
		return geneName;
	}

	public void setGeneName(List<String> geneName) {
		this.geneName = geneName;
	}

	public double getpValue() {
		return pValue;
	}

	public void setpValue(double pValue) {
		this.pValue = pValue;
	}
}

