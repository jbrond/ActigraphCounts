package timo.jyu;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import timo.jyu.filter.ButterworthCoefficients;

public class Utils{
	//Returns 4th order band pass butterworth filter coefficients
	//http://www-users.cs.york.ac.uk/~fisher/mkfilter/
	public static double[][] getBandPassButterworthCoefficients(double[] corners, double sRate){
		ButterworthCoefficients b = new ButterworthCoefficients();
		b.butter(new String[]{"Bu","Bp","o","4","a",Double.toString(0.1/sRate),Double.toString(7d/sRate)});
		double[][] ret = new double[2][];
		ret[0] = new double[b.xcoeffs.length];
		ret[1] = new double[b.ycoeffs.length];
		//Divide the xcoeffs (b in matlab) with gain, and reverse the order of ycoeffs(a in matlab)
		for (int i = 0;i<b.xcoeffs.length;++i){
			ret[0][i] =b.xcoeffs[i]/b.gain;
			ret[1][b.ycoeffs.length-1-i] =-1d*b.ycoeffs[i];
		}
		return ret;
	}
	
	private static double[] reverse(double[] a){
		double[] b = new double[a.length];
		for (int i = 0;i<b.length;++i){
			b[b.length-1-i] = a[i];
		}
		return b;		
	}
	
	private static double sum(double[] a){
		double b = 0d;
		for (int i = 0;i<a.length;++i){
			b += a[i];
		}
		return b;
	}
	/**helper function for zero-lag filtering
		copied from octave filtfilt.m*/
	public static double[] prepState(double[] b, double[] a,double[] state){
		double kdc = sum(b)/sum(a);
		double[] temp = new double[a.length];
		//Initialise state if kdc is not infinity or NaN
		if (Math.abs(kdc) < Double.POSITIVE_INFINITY && !Double.isNaN(kdc)){
			
			for (int i = 0;i<a.length;++i){
				temp[i] = b[i]-kdc*a[i];
			}
			temp = reverse(temp);
			for (int i = 1;i<temp.length;++i){
				temp[i] += temp[i-1];
			}
			temp = reverse(temp);
			for (int i = 1;i<temp.length;++i){
				state[i-1] = temp[i];
			}
		}
		return state;
	}
	
	private static double[] multArray(double[] a, double b){
		double[] c = new double[a.length];
		for (int i = 0;i<a.length;++i){c[i] = b*a[i];}
		return c;
	}
	
	public static double[] getMirrored(double[] signal, int mirrorLength){
		double[] temp = new double[signal.length+2*mirrorLength];
		for (int i = 0; i<mirrorLength; ++i){
			temp[i] = 2d*signal[0]-signal[mirrorLength-i];
			temp[mirrorLength+signal.length+i]=2d*signal[signal.length-1]-signal[signal.length-2-i];
		}
		for (int i = 0;i<signal.length;++i){
			temp[mirrorLength+i] = signal[i];
		}
		return temp;
	}
	
	/**
		Zero-lag filter. Matches Octave/Matlab filtfilt 
	*/
	public static double[] filtfilt(double[] b, double[] a, double[] signal){
		double[] state = prepState(b,a,new double[a.length-1]);
		//Mirror signal from the beginning and from the end, and insert signal in the middle
		int initBackwardSamples = signal.length < 3*(a.length-1) ? signal.length : 3*(a.length-1);
		double[] temp = getMirrored(signal,initBackwardSamples);
		temp = filter(b,a,temp,multArray(state,temp[0]));	//Filter forward
		temp = reverse(temp);	//Switch direction
		temp = filter(b,a,temp,multArray(state,temp[0]));	//Filter backward
		temp = reverse(temp);	//Switch direction
		//Return the mid-part without mirrored data
		double[] output = new double[signal.length];
		for (int i = 0;i<signal.length;++i){
			output[i] = temp[initBackwardSamples+i];
		}
		return output;
	}
	
	
	/**
		Filter with initial state
		ported from https://github.com/greenm01/forc/blob/master/sandbox/filter.m
	*/
	public static double[] filter(double[] b, double[] a, double[] signal,double[] state){
		double[] output = new double[signal.length];
		for (int i = 0;i<output.length;++i){
			output[i] = state[0]+b[0]*signal[i]; //Next output does not depend on a coeffs
			//Update state vector after the current output is known
			for (int j = 1; j<state.length;++j){
				state[j-1] = state[j]-a[j]*output[i]+b[j]*signal[i];
			}
			state[state.length-1]=b[state.length]*signal[i]-a[state.length]*output[i];		
		}
		return output;
	}
	
		/**
		Filter without initial state (=set initial state to all zeros)
	*/
	public static double[] filter(double[] b, double[] a, double[] signal){
		return filter(b,a,signal,new double[a.length-1]);
	}

	public static double[] interp(double[] x, double[] y, double[] xx){
		double[] yy = new double[xx.length];
		PolynomialSplineFunction interpolator = (new LinearInterpolator()).interpolate(x, y);
		for (int i = 0;i<yy.length;++i){
			yy[i] = interpolator.value(xx[i]);
		}
		return yy;
	}
	
}
