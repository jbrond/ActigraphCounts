/*
	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	N.B.  the above text was copied from http://www.gnu.org/licenses/gpl.html
	unmodified. I have not attached a copy of the GNU license to the source...

    Copyright (C) 2011 - 2019 Timo Rantalainen, tjrantal@gmail.com
*/

package timo.jyu.filter;
/*Modified by Timo Rantalainen 2011 to 2019 from http://www-users.cs.york.ac.uk/~fisher/mkfilter/*/
/* mkfilter -- given n, compute recurrence relation
   to implement Butterworth, Bessel or Chebyshev filter of order n
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   September 1992 */
public class ButterworthCoefficients{

	public static final int opt_bu =  0x00002;	/* -Bu		Butterworth characteristic     */
	public static final int opt_lp = 0x00020;	/* -Lp		lowpass			       */
	public static final int opt_hp = 0x00040;	/* -Hp		highpass		       */
	public static final int opt_bp = 0x00080;	/* -Bp		bandpass		       */
	public static final int opt_a = 0x00400;	/* -a		alpha value		       */
	public static final int opt_o = 0x01000;	/* -o		order of filter		       */

	public class Pzrep{ 		
		public Complex[] poles;
		public Complex[] zeros;
		public int numpoles;
		public int numzeros;
		public Pzrep(int maxpz){
			 poles = new Complex[maxpz];
			 zeros = new Complex[maxpz];
			 numpoles = 0;
			 numzeros = 0;
		}
		
	};

	public Pzrep splane;
	public Pzrep zplane;
	public int order;
	public double raw_alpha1, raw_alpha2, raw_alphaz;
	public Complex dc_gain, fc_gain, hf_gain;
	public int options;
	public double warped_alpha1, warped_alpha2, chebrip, qfactor;
	public int polemask;
	public double[] xcoeffs, ycoeffs, xcoeff,ycoeff;
	public double gain;
	public int nx, ny;
	int MAXPZ;
	
	public ButterworthCoefficients(){
		
	}
	
	public static void main(String[] args){
		 ButterworthCoefficients butterworthCoefficients = new  ButterworthCoefficients();
		 butterworthCoefficients.butter(args);
	}
	
	public void butter(String[] argv)
	  { //System.out.println("Script started");

		readcmdline(argv);
		//System.out.println("Command line read "+options+ " o "+order+" lp "+raw_alpha1);
		MAXPZ = 2*order;
		splane = new Pzrep(MAXPZ);
		zplane = new Pzrep(MAXPZ);
		xcoeffs = new double[MAXPZ+1];
		ycoeffs = new double[MAXPZ+1];

		
		//checkoptions();
		//printf("Saatiin optiot\n");
		setdefaults();
		//System.out.println("Defaults set");
		
		compute_s();
	    prewarp();
	    normalize();
		compute_z_blt();
		
		//printf("Polynomia kasvattamaan\n");
		//System.out.println("ZNumzeros "+ zplane.numzeros+" ZNumPoles "+zplane.numpoles);
		//System.out.println("raw alpha "+ raw_alpha1+" warped alpha "+warped_alpha1);
		expandpoly();
		//System.out.println("Polynome expanded");
		Complex gains = new Complex(0.0);
		//System.out.println("DC re "+dc_gain.re+" im" + dc_gain.im);
		if ((options & opt_lp)> 0) gains = dc_gain;
		if ((options & opt_hp)> 0) gains = hf_gain;
		if ((options & opt_bp)> 0) gains = fc_gain;
	   //printf("Saatiin gainit\n");
		gain = Complex.hypot(gains);
		////System.out.println("Gain "+gain);
		nx = zplane.numzeros;
		ny = zplane.numpoles;
		//printresults(argv);

	  }

	  
	  /*Additional functions*/
	public void prewarp()
	  { /* for bilinear transform, perform pre-warp on alpha values */
		warped_alpha1 = Math.tan(Math.PI * raw_alpha1) / Math.PI;
		warped_alpha2 = Math.tan(Math.PI * raw_alpha2) / Math.PI;
	  }
	  /*Additional functions end*/
	  
	  
	  
	public void readcmdline(String[] argv)
	  { options = order = polemask = 0;
		int ap = 0;
		//printf("Mennaan katsomaan montako on optioita\n");
		while (ap < argv.length)
		  { 
			  //printf("Optio No. %d\n",ap);
			  int m = decodeoptions(argv[ap++]);
			if ((m & opt_a)>0){ 
				raw_alpha1 = getfarg(argv[ap++]);
				 
				if (ap < argv.length && argv[ap].charAt(0) != '-') {
					raw_alpha2 = getfarg(argv[ap++]);
				} else {
					raw_alpha2 = raw_alpha1;
				}
			}
			if ((m & opt_o)>0) order = getiarg(argv[ap++]);
			options |= m;
		  }
	  }

	public boolean seq(String s1, String s2){
		return s1.compareTo(s2) == 0;
	}

	public int decodeoptions(String s)
	  { 
		int m = 0;
		if (seq(s,new String("Bu"))) m |= opt_bu;
		else if (seq(s,new String( "Lp"))) m |= opt_lp;
		else if (seq(s,new String( "Hp"))) m |= opt_hp;
		else if (seq(s,new String( "Bp"))) m |= opt_bp;

		else
		  { 
			  int ind = 0;
			  int keepOn = 1;
			  while (keepOn != 0)
			  { int bit = optbit(s.charAt(ind));
				if (bit == 0){ 
					//System.out.println("Erroneous commandline options");				
				}
				m |= bit;
				if (bit != 0 || ind >= s.length()){
					keepOn = 0;
				}
			  }
		  }
		return m;
	  }

	static int optbit(char c)
	  { switch (c)
		  { 
			case 'a':   return opt_a;
			case 'o':   return opt_o;
			default:    return 0;	
		  }
	  }

	public double getfarg(String s)
	  { if (s.isEmpty()) {
			//System.out.println("Error in getfarging");
		}
		return Double.parseDouble(s);
	  }

	public int getiarg(String s)
	  { if (s.isEmpty()){
		//System.out.println("Error in getiarging");
		}
		return Integer.parseInt(s);
	  }


	public void setdefaults()
	  { 
		polemask = 0xFFFFFFFF;
		if ((options & opt_bp) == 0) raw_alpha2 = raw_alpha1;
	  }

	public void compute_s() /* compute S-plane poles for prototype LP filter */
	  { 
		////System.out.println("Compute s");
		
		splane.numpoles = 0;
		////System.out.println("If...");
		if ((options & opt_bu)> 0)
		  { /* Butterworth filter */
			////System.out.println("Butterworth filter");
			for (int i = 0; i < 2*order; i++)
			  { double theta = 0;
				if ((order & 1) > 0) {
					////System.out.println("Pariton");
					//theta = (((double) i)*Math.PI) / (double) order;
					}
				if ((order & 1) == 0){
					
					theta = ((((double) i)+0.5)*Math.PI) / (double) order;
					////System.out.println("Parillinen "+theta );
				}
				choosepole(Complex.expj(theta));
				////System.out.println("polemask "+polemask+" theta "+ theta);
			  }
		  }
	  }

	public void choosepole(Complex z){ 
		
		if (z.re < 0.0){ 
			////System.out.println("re "+z.re+" im "+z.im);
			if ((polemask & 1) > 0){ 
				splane.poles[splane.numpoles++] = z;
				////System.out.println("PoleZhosen "+splane.numpoles+" re "+z.re+" im "+z.im);
			
			}
			polemask = polemask >> 1;
		}
	}

	public void normalize()		/* called for trad, not for -Re or -Pi */
	  { double w1 = 2*Math.PI * warped_alpha1;
		double w2 = 2*Math.PI * warped_alpha2;
		/* transform prototype into appropriate filter type (lp/hp/bp/bs) */
		switch (options & (opt_lp | opt_hp | opt_bp)){ 
			case opt_lp:
			  { for (int i = 0; i < splane.numpoles; i++){ splane.poles[i] = Complex.multiplyWithConstant(w1,splane.poles[i]);
				////System.out.println("Warped s-poles re "+splane.poles[i].re+" im "+splane.poles[i].im);
			  }
				splane.numzeros = 0;
				
				break;
			  }

			case opt_hp:
			  { int i;
				for (i=0; i < splane.numpoles; i++) splane.poles[i] = Complex.divide(new Complex(w1),splane.poles[i]);
				for (i=0; i < splane.numpoles; i++) splane.zeros[i] = new Complex(0.0);	 /* also N zeros at (0,0) */
				splane.numzeros = splane.numpoles;
				break;
			  }

			case opt_bp:
			  { double w0 = Math.sqrt(w1*w2), bw = w2-w1; int i;
				for (i=0; i < splane.numpoles; i++)
				  { Complex hba = Complex.multiplyWithConstant(0.5,Complex.multiplyWithConstant(bw,splane.poles[i]));
				Complex temp = Complex.csqrt(Complex.subtract(new Complex(1.0),Complex.multiply(Complex.divide(new Complex(w0),hba),Complex.divide(new Complex(w0),hba))));
				splane.poles[i] = Complex.multiply(hba,Complex.sum(new Complex(1.0),temp));
				splane.poles[splane.numpoles+i] = Complex.multiply(hba,(Complex.subtract(new Complex(1.0),temp)));
				  }
				for (i=0; i < splane.numpoles; i++) splane.zeros[i] = new Complex(0.0);	 /* also N zeros at (0,0) */
				splane.numzeros = splane.numpoles;
				splane.numpoles *= 2;
				break;
			  }

		  }
	  }

	public void compute_z_blt() /* given S-plane poles & zeros, compute Z-plane poles & zeros, by bilinear transform */
	  { int i;
		zplane.numpoles = splane.numpoles;
		zplane.numzeros = splane.numzeros;
		for (i=0; i < zplane.numpoles; i++){
			zplane.poles[i] = blt(splane.poles[i]);
			////System.out.println("Warped z-poles re "+zplane.poles[i].re+" im "+zplane.poles[i].im);
		}
		for (i=0; i < zplane.numzeros; i++){ zplane.zeros[i] = blt(splane.zeros[i]);}
		while (zplane.numzeros < zplane.numpoles) zplane.zeros[zplane.numzeros++] = new Complex(-1.0);
	  }

	public Complex blt(Complex pz)
	  { return Complex.divide(Complex.sum(new Complex(2.0),pz),Complex.subtract(new Complex(2.0),pz));
	  }

	public void compute_z_mzt() /* given S-plane poles & zeros, compute Z-plane poles & zeros, by matched z-transform */
	  { int i;
		zplane.numpoles = splane.numpoles;
		zplane.numzeros = splane.numzeros;
		for (i=0; i < zplane.numpoles; i++) zplane.poles[i] = Complex.cexp(splane.poles[i]);
		for (i=0; i < zplane.numzeros; i++) zplane.zeros[i] = Complex.cexp(splane.zeros[i]);
	  }




	public void expandpoly() /* given Z-plane poles & zeros, compute top & bot polynomials in Z, and then recurrence relation */
	  { Complex[] topcoeffs = new Complex[MAXPZ+1];
		Complex[] botcoeffs = new Complex[MAXPZ+1];
		int i;
		expand(zplane.zeros, zplane.numzeros, topcoeffs);
		expand(zplane.poles, zplane.numpoles, botcoeffs);
		dc_gain = Complex.evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, new Complex(1.0));
		//System.out.println("dc_gain "+dc_gain.re+" "+dc_gain.im);
		double theta = 2.0*Math.PI * 0.5 * (raw_alpha1 + raw_alpha2); /* "jwT" for centre freq. */
		fc_gain = Complex.evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, Complex.expj(theta));
		hf_gain = Complex.evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, new Complex(-1.0));
		for (i = 0; i <= zplane.numzeros; i++) {
			xcoeffs[i] = +(topcoeffs[i].re / botcoeffs[zplane.numpoles].re);
			//System.out.println("xcoeff "+i+" "+xcoeffs[i]);
		}
		for (i = 0; i <= zplane.numpoles; i++) {
			ycoeffs[i] = -(botcoeffs[i].re / botcoeffs[zplane.numpoles].re);
			//System.out.println("ycoeff "+i+" "+ycoeffs[i]);
		}
	  }

	public void expand(Complex[] pz, int npz, Complex[] coeffs)
	  { /* compute product of poles or zeros as a polynomial of z */
		int i;
		coeffs[0] = new Complex(1.0);
		for (i=0; i < npz; i++) coeffs[i+1] = new Complex(0.0);
		for (i=0; i < npz; i++) multin(pz[i], npz, coeffs);
	  }

	public void multin(Complex w, int npz, Complex[] coeffs)
	  { /* multiply factor (z-w) into coeffs */
		Complex nw = Complex.changeSign(w);
		////System.out.println("W "+w.re+" "+w.im);
		////System.out.println("-W "+nw.re+" "+nw.im);
		for (int i = npz; i >= 1; i--) coeffs[i] = Complex.sum(Complex.multiply(nw,coeffs[i]),coeffs[i-1]);
		coeffs[0] = Complex.multiply(nw,coeffs[0]);
	  }

	public void printresults(String[] argv)
	  { 
		Complex gain = new Complex(0.0);
		if ((options & opt_lp) > 0) gain = dc_gain;
		if ((options & opt_hp) > 0) gain = hf_gain;
		if ((options & opt_bp) > 0) gain = fc_gain;

		System.out.println("Gain  = "+ Complex.hypot(gain)+" "+ Complex.hypot(gain));
		printcoeffs(new String("NZ"), zplane.numzeros, xcoeffs);
		printcoeffs(new String("NP"), zplane.numpoles, ycoeffs);

	  }


	public void printcoeffs(String pz, int npz, double coeffs[])
	  { //System.out.println(pz +" N = "+ npz);
		for (int i = 0; i <= npz; i++){  
			System.out.println(coeffs[i]);
		}
	  }

	public double[] filtfilt(double[] data){
		//System.out.println("nx "+nx+" ny "+ny+" xcoeff "+xcoeff.length+" ycoeff "+ycoeff.length);
		int dataLength = data.length;
		//double[] filtered = new double[dataLength];
		double[] xv = new double[nx+1];
		double[] yv = new double[ny+1];
		//Täytellään xv ja yv ensimmäisellä arvolla
		for (int j = 0;j<=nx;j++){
			xv[j]= data[0]/gain;
		}

		for (int j = 0;j<=ny;j++){
			yv[j]= data[0];
		}
		//Filtteröidään etuperin
		for (int i = 0;i<dataLength;i++){
			for (int j = 0;j<nx;j++){
				xv[j]= xv[j+1];
			}
			xv[nx] = data[i]/gain;
			for (int j = 0;j<ny;j++){
				yv[j]= yv[j+1];
			}
			yv[ny] = 0;
			for (int j = 0; j<=nx;j++){
				yv[ny] = yv[ny]+xv[j]*xcoeffs[j];
			}
			
			for (int j = 0;j<ny;j++){
				yv[ny]=yv[ny]+yv[j]*ycoeffs[j];
			}
			data[i] = yv[ny];
		}
		//Filtteröinti takaperin
		//Täytellään xv ja yv viimeisellä arvolla
		for (int j = 0;j<=nx;j++){
			xv[j]= data[dataLength-1]/gain;
		}

		for (int j = 0;j<=ny;j++){
			yv[j]= data[dataLength-1];
		}
		//Filtterointi takaperin
		for (int i = dataLength-1;i>=0;i--){
			for (int j = 0;j<nx;j++){
				xv[j]= xv[j+1];
			}
			xv[nx] = data[i]/gain;
			for (int j = 0;j<ny;j++){
				yv[j]= yv[j+1];
			}
			yv[ny] = 0;
			for (int j = 0; j<=nx;j++){
				yv[ny] = yv[ny]+xv[j]*xcoeffs[j];
			}
			for (int j = 0;j<ny;j++){
				yv[ny]=yv[ny]+yv[j]*ycoeffs[j];
			}
			data[i] = yv[ny];
		}
		return data;
	}

}
