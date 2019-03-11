import java.util.List;

public class MergingPojo {
	List<String> mergedPath;
	int pathNumber_2;
	int pathNumber_1;
	Double combineCP;

	public Boolean compareMergingPojosObject(MergingPojo o)
	{
		Boolean still = false;
		if (pathNumber_2==o.pathNumber_2)
			still = true;
		else
			still = false;
		if (still)
			if (pathNumber_1==o.pathNumber_1)
				still = true;
			else
				still = false;
		if (still)
		{
			if (mergedPath.equals(o.mergedPath))
				still = true;
			else
					still = false;

		}
		return still;
	}

	public int getPathNumber_2() {
		return pathNumber_2;
	}
	public void setPathNumber_2(int pathNumber_2) {
		this.pathNumber_2 = pathNumber_2;
	}
	public int getPathNumber_1() {
		return pathNumber_1;
	}
	public void setPathNumber_1(int pathNumber_1) {
		this.pathNumber_1 = pathNumber_1;
	}
	public List<String> getMergedPath() {
		return mergedPath;
	}
	public void setMergedPath(List<String> mergedPath) {
		this.mergedPath = mergedPath;
	}
	public Double getCombineCP() {
		return combineCP;
	}
	public void setCombineCP(Double combineCP) {
		this.combineCP = combineCP;
	}

}
