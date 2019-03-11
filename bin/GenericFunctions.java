
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import java.math.BigDecimal;

public class GenericFunctions {

	static double hartungFunction(List<Double> pvalueList) {
		if (pvalueList.size() == 1 || pvalueList == null)
			return Double.NaN;
		List<BigDecimal> qnormList = new ArrayList<>();
		int n = pvalueList.size();
			for (Double pavlue : pvalueList) {
				try{
				qnormList.add(new BigDecimal(NormalCDFInverse(pavlue)));
				}catch (Exception e)
				{
					return Double.NaN;
				}
			}
		BigDecimal sum = new BigDecimal(0.0);
		BigDecimal avt = new BigDecimal(0.0);
		for (BigDecimal qnorm : qnormList) {
			sum = sum.add(qnorm);
		}
		BigDecimal num = new BigDecimal(n);
		if (num.compareTo(BigDecimal.ZERO) == 0) {
			return Double.NaN;
		}
		avt = sum.divide(num, 5, RoundingMode.HALF_UP);
		BigDecimal q = new BigDecimal(0.0);
		BigDecimal loopSum = new BigDecimal(0.0);
		BigDecimal minusOneBD = new BigDecimal(-1);
		for (BigDecimal double1 : qnormList) {
			BigDecimal tempValue = double1.subtract(avt);
			if (tempValue.compareTo(BigDecimal.ZERO) == -1)
				tempValue = tempValue.multiply(minusOneBD);
			BigDecimal powerResult = customPowerFucntion(tempValue, 2);
			loopSum = loopSum.add(powerResult);
		}
		int xx = n - 1;
		q = loopSum.divide(BigDecimal.valueOf(xx), 5, RoundingMode.HALF_UP);
		BigDecimal one = new BigDecimal(1.0);
		BigDecimal rhohat = one.subtract(q);
		BigDecimal rhostar = new BigDecimal(0.0);
		BigDecimal nDouble = new BigDecimal(n);
		BigDecimal divideRresult = new BigDecimal(n);
		BigDecimal minus1 = new BigDecimal(-1.0);
		BigDecimal subtractResult = new BigDecimal(n);
		subtractResult = nDouble.subtract(one);
		divideRresult = one.divide(subtractResult, 5, RoundingMode.HALF_UP);
		divideRresult = divideRresult.multiply(minus1);
		rhostar = maxfunction(divideRresult, rhohat);
		double kappa = 0.0;
		kappa = .2;
		List<Integer> lambdaList = new ArrayList<Integer>();
		for (Double pvalue : pvalueList) {
			lambdaList.add(1);
		}
		double tempAnser = primitiveFunction(n, qnormList, rhostar.doubleValue(), kappa, lambdaList);
		return tempAnser;
	}

	static BigDecimal maxfunction(BigDecimal a, BigDecimal b) {
		if (a.compareTo(b) == 1)
			return a;
		return b;
	}

	static double calcualteKappa(int n, double rostat) {
		int count = 1 + 1 / (n - 1);
		return (count - rostat) / 10;
	}

	static double primitiveFunction(int n, List<BigDecimal> tvalues, double rhostar, double kappa,
			List<Integer> lambdaList) {
		kappa = .2;
		BigDecimal sumofLambdaTValyes = new BigDecimal(0.0);
		BigDecimal LambdaPowerTwoSum = new BigDecimal(0.0);
		List<BigDecimal> lambdaPowerTwoList = new ArrayList<>();
		List<BigDecimal> multipleLamdatValues = new ArrayList<>();
		BigDecimal twoBD = new BigDecimal(2);
		for (int xx = 0; xx < tvalues.size(); xx++) {
			BigDecimal temp = tvalues.get(xx).multiply(new BigDecimal(lambdaList.get(xx)));
			BigDecimal lambdaPower = BigDecimalMath.pow(new BigDecimal(lambdaList.get(xx)), twoBD);
			lambdaPowerTwoList.add((lambdaPower));
			multipleLamdatValues.add(temp);
		}
		for (BigDecimal double1 : multipleLamdatValues) {
			BigDecimal tempX = new BigDecimal(0.0);
			tempX = sumofLambdaTValyes.add(double1);
			sumofLambdaTValyes = tempX;
		}
		for (BigDecimal lasmdapower : lambdaPowerTwoList) {
			BigDecimal tempX = new BigDecimal(0.0);
			tempX = LambdaPowerTwoSum.add(lasmdapower);
			LambdaPowerTwoSum = tempX;
		}
		BigDecimal partOneAnswer = sumofLambdaTValyes;
		BigDecimal PowerTwoOfLambdaSum = new BigDecimal(0.0);
		int tempAddition = 0;
		for (Integer lambda : lambdaList) {
			tempAddition += lambda;
		}
		PowerTwoOfLambdaSum = new BigDecimal(Math.pow(tempAddition, 2));
		BigDecimal partTwoAnser = PowerTwoOfLambdaSum.subtract(LambdaPowerTwoSum);
		double doubleNN = n;
		BigDecimal srqt = new BigDecimal(Math.sqrt(2 / (doubleNN - 1)));
		BigDecimal OneMinusRhoStar = new BigDecimal(1 - rhostar);
		BigDecimal kappaBD = new BigDecimal(kappa);
		BigDecimal kapasrqt = kappaBD.multiply(srqt);
		BigDecimal mult1 = kapasrqt.multiply(OneMinusRhoStar);
		BigDecimal partThreeAnser = mult1.add(new BigDecimal(rhostar));
		BigDecimal p2p3 = partTwoAnser.multiply(partThreeAnser);
		BigDecimal Lambdap2p3 = LambdaPowerTwoSum.add(p2p3);
		Double xxMatg = Math.sqrt(Lambdap2p3.doubleValue());
		BigDecimal combineAner = partOneAnswer.divide(new BigDecimal(xxMatg), 3, RoundingMode.HALF_UP);
		BigDecimal finalAnserBD = pnormBD2(combineAner);
		return finalAnserBD.doubleValue();
	}

	static double repFunction(int length) {
		return 0.0;
	}

	static double RationalApproximation(double t) {

		double c[] = { 2.515517, 0.802853, 0.010328 };
		double d[] = { 1.432788, 0.189269, 0.001308 };
		return t - ((c[2] * t + c[1]) * t + c[0]) / (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
	}

	static double NormalCDFInverse(double p) {
		if (p <= 0.0 || p >= 1.0) {
              return 0;
		}
		if (p < 0.5) {
			return -RationalApproximation(Math.sqrt(-2.0 * Math.log(p)));
		} else {
			return RationalApproximation(Math.sqrt(-2.0 * Math.log(1 - p)));
		}
	}

	static BigDecimal customPowerFucntion(BigDecimal a, int b) {
		BigDecimal minusOne = new BigDecimal(-1);
		if (a.compareTo(minusOne) == 0)
			a = a.multiply(minusOne);
		BigDecimal result = new BigDecimal(1);
		for (int i = 1; i <= b; i++) {
			result = result.multiply(a);
		}
		DecimalFormat df = new DecimalFormat("0.000000");
		BigDecimal returnResult = new BigDecimal(0);
		String formateNormal = df.format(result);
		try {
			returnResult = (new BigDecimal((double) df.parse(formateNormal)));
		} catch (Exception e) {
		}
		return returnResult;
	}

	static BigDecimal pnormBD2(BigDecimal xBD) {

		BigDecimal zero = new BigDecimal(0.0);
		BigDecimal one = new BigDecimal(1.0);
		BigDecimal neg = new BigDecimal(0.0);
		BigDecimal negOne = new BigDecimal(-1.0);
		if (zero.compareTo(xBD) == 1)
			neg = new BigDecimal(1.0);
		if (neg.compareTo(one) == 0) {
			xBD = xBD.multiply(negOne);
		}
		BigDecimal point2X = xBD.multiply(new BigDecimal(0.2316419));
		BigDecimal point2XPlus1 = one.add(point2X);
		BigDecimal k = one.divide(point2XPlus1, 5, RoundingMode.UP);
		BigDecimal onepoint3wihk = k.multiply(new BigDecimal(1.330274429));
		BigDecimal onepoint3wih3Minus1pint8 = onepoint3wihk.subtract(new BigDecimal(1.821255978));
		BigDecimal onepoint3wih3Minus1pint8KK = k.multiply(onepoint3wih3Minus1pint8);
		BigDecimal onepoint3wih3Minus1pint8KKOnePint7 = onepoint3wih3Minus1pint8KK.add(new BigDecimal(1.781477937));
		BigDecimal pointOnePint7KK = onepoint3wih3Minus1pint8KKOnePint7.multiply(k);
		BigDecimal pointOnePint7KKMinusP3 = pointOnePint7KK.subtract(new BigDecimal(0.356563782));
		BigDecimal pointOnePint7KKMinusP3KK = pointOnePint7KKMinusP3.multiply(k);
		BigDecimal pointOnePint7KKMinusP3KKAddPoint3 = pointOnePint7KKMinusP3KK.add(new BigDecimal(0.319381530));
		BigDecimal y = pointOnePint7KKMinusP3KKAddPoint3.multiply(k);
		double minuspoint5 = (-.5);
		double xbd = xBD.doubleValue();
		double xMulX = xbd * xbd;
		double xMulXPoint5 = xMulX * (minuspoint5);
		double Euler = 2.7182818;
		double minusOneBD = (-1);
		if (xMulXPoint5 == -1)
			xMulXPoint5 = xMulXPoint5 * (minusOneBD);
		double mathEx = Math.pow(Euler, xMulXPoint5);
		BigDecimal Ex = BigDecimal.valueOf(mathEx);
		BigDecimal mathExMulY = Ex.multiply(y);
		BigDecimal mathExMulYMulPoint3 = mathExMulY.multiply(new BigDecimal(0.398942280401));
		BigDecimal yF = one.subtract(mathExMulYMulPoint3);
		BigDecimal oneMinusNeg = one.subtract(neg);
		BigDecimal oneMinusY = one.subtract(yF);
		BigDecimal megMultioneMinusY = neg.multiply(oneMinusY);
		BigDecimal OneMinusNeg = oneMinusNeg.multiply(yF);
		BigDecimal finalAddition = OneMinusNeg.add(megMultioneMinusY);

		return finalAddition;
	}
}
