package timo.jyu;

/*
	Java port of Jan Brønd's Actigraph counts from raw accelerations
	Written by Timo Rantalainen tjrantal at gmail dot com 2019
	This Java port released to the public domain with a CC0-BY Creative Commons Attribution license	
	Filtering ported from http://www-users.cs.york.ac.uk/~fisher/mkfilter/, https://github.com/greenm01/forc/blob/master/sandbox/filter.m and https://searchcode.com/codesearch/view/9521190/ 
*/

public class ActigraphCounts{
	//Coeffs written to a text file from Matlab
	private static final double[] A = new double[]{1,-4.16372602554363,7.57115309014007,-7.9804690250911,5.38501191026769,-2.46356271321257,0.89238142271725,0.0636099868633388,-1.34810512714076,2.47338133053049,-2.92571735841718,2.92983230386598,-2.78159062882719,2.4776735357121,-1.68473849390463,0.464828627239016,0.465652889035618,-0.67311896742996,0.416203225759379,-0.13832322391961,0.0198517159761605};
	private double[] B = new double[]{0.0491089825140488,-0.122841835307157,0.143557884896153,-0.112693989220238,0.0538037410952924,-0.020230273840001,0.00637784647673757,0.0185125409235852,-0.0381541058906574,0.0487265187117185,-0.0525772146919335,0.047847138089546,-0.0460148280299714,0.0362833364868511,-0.0129768121654561,-0.00462621079355692,0.0128353970741233,-0.00937622141658307,0.00344850106651387,-0.000809720155277696,-0.000196225290878896};
	public static final double adcResolution = 0.0164;
	public static final double gain = 0.965;
	public static final double deadband = 0.068;
	public static final double maxTrunc = 2.13;
	public static final int integN = 10;
	double[] counts = null;
	double[][] afCoeffs = null;
		
	/*Constructor that runs the analysis*/
	public ActigraphCounts(double[] resultant, double sRate){
		afCoeffs =  Utils.getBandPassButterworthCoefficients(new double[]{0.1,7d}, sRate);
		for (int i=0;i<B.length;++i){
			B[i] = B[i]*gain; //Apply ag filter gain on B coeffs
		}	
		counts = calcCounts(pptruncDeadband8bit(down10(getAGFiltered(getInterpolated(getAFiltered(resultant),sRate))),maxTrunc,deadband,adcResolution));
	}
	
	/*Constructor to be used to run the analysis step by step*/
	public ActigraphCounts(double sRate){
		afCoeffs =  Utils.getBandPassButterworthCoefficients(new double[]{0.1,7d}, sRate);
		for (int i=0;i<B.length;++i){
			B[i] = B[i]*gain; //Apply ag filter gain on B coeffs
		}	
	}
	
	public double[] getCounts(double[] a){
		counts = calcCounts(a);
		return counts;
	}
	
	public void setA(double[] a){
		afCoeffs[1] = new double[a.length];
		for (int i = 0; i<a.length; ++i){afCoeffs[1][i] = a[i];}
	}
	
	public void setB(double[] a){
		afCoeffs[0] = new double[a.length];
		for (int i = 0; i<a.length; ++i){afCoeffs[0][i] = a[i];}
	}
	
	public double[] getA(){
		return afCoeffs[1];
	}
	
	public double[] getB(){
		return afCoeffs[0];
	}
	
	public double[] getCounts(){
		return counts;
	}
	
	private double[] getT(int length, double sRate){
		double[] t = new double[length];
		for (int i = 0; i< length; ++i){
			t[i] = ((double) i)/sRate;
		}
		return t;
	}
	
	public double[] getInterpolated(double[] afiltered,double sRate){
		//Resample the signal to 30 Hz
		double[] t = getT(afiltered.length,sRate);
		double[] intT = getT((int) (Math.floor(t[t.length-1]*30d))+1,30d);
		return Utils.interp(t,afiltered,intT);
	}
	
	public double[] getAFiltered(double[] resultant){
		return Utils.filtfilt(afCoeffs[0],afCoeffs[1],resultant);
	}
	
	public double[] getAGFiltered(double[] down30){
		return Utils.filter(B,A,down30);
	}
	
	//AG counts are calculated based on 10 Hz data. Downsample 30 Hz to 10 Hz
	//By taking every 3rd sample
	public double[] down10(double[] a){
		double[] b = new double[(int) (Math.floor(((double)a.length)/3d))];
		for (int i = 0; i< b.length; ++i){
			b[i] = a[3*i];
		}
		return b;
	}
	
	public double[] pptruncDeadband8bit(double[] a, double maxMagnitude, double dbMagnitude,double adb){
		double[] b = new double[a.length];
		for (int i = 0; i< b.length; ++i){
			b[i] = Math.abs(a[i]) <=maxMagnitude ? a[i] : Math.signum(a[i])*maxMagnitude;	//Truncate max val
			b[i] = Math.abs(b[i]) >=dbMagnitude ? Math.abs(b[i]) : 0d;	//Apply deadband and rectify
			b[i] = Math.floor(b[i]/adb);	//8-bit values
		}
		return b;
	}
	
	//Calculate 10 sample integrals
	public double[] calcCounts(double[] a){
		double[] b = new double[(int) (Math.ceil(((double) a.length)/((double) integN)))];
		for (int i = 0;i<b.length;++i){
			for (int j = 0;j<integN;++j){
				if (integN*i+j <a.length){
					b[i]+=a[integN*i+j];
				}
			}
		}
		return b;	
	}
}
