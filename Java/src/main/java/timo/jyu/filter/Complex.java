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

    Copyright (C) 2011- 2019 Timo Rantalainen, tjrantal@gmail.com
*/

package timo.jyu.filter;

/*Modified by Timo Rantalainen 2011 to 2019 from http://www-users.cs.york.ac.uk/~fisher/mkfilter/*/

/* mkfilter -- given n, compute recurrence relation
   to implement Butterworth, Bessel or Chebyshev filter of order n
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   September 1992 */

/* Routines for complex arithmetic */
public class Complex{
	public double re;
	public double im;
	/*Constructors*/
	public Complex(double reIn, double imIn){
		re = reIn;
		im = imIn;
	}
	
	public Complex(double reIn){
		re = reIn;
		im = 0;
	}
	/*Constructors done*/
	
	/*methods*/
	public static Complex evaluate(Complex[] topco, int nz, Complex[] botco, int np, Complex z)
	  { /* evaluate response, substituting for z */
		Complex divident = eval(topco, nz, z);
		Complex divisor = eval(botco, np, z);
		return divide(divident,divisor);
	  }

	public static Complex eval(Complex[] coeffs, int npz, Complex z)
	  { /* evaluate polynomial in z, substituting for z */
		Complex result = new Complex(0.0);
		for (int i = npz; i >= 0; i--) result = sum((multiply(result,z)),coeffs[i]);
		return result;
	  }

	public static double hypot(Complex z) { return Math.hypot(z.im, z.re); }
	public static double atan2(Complex z) { return Math.atan2(z.im, z.re); }
	  
	public static Complex csqrt(Complex x)
	  { double r = hypot(x);
		Complex z = new Complex(Xsqrt(0.5 * (r + x.re)),
				Xsqrt(0.5 * (r - x.re)));
		if (x.im < 0.0) z.im = -z.im;
		return z;
	  }

	public static double Xsqrt(double x)
	  { /* because of deficiencies in hypot on Sparc, it's possible for arg of Xsqrt to be small and -ve,
		   which logically it can't be (since r >= |x.re|).	 Take it as 0. */
		return (x >= 0.0) ? Math.sqrt(x) : 0.0;
	  }

	public static Complex cexp(Complex z)
	  { return multiplyWithConstant(Math.exp(z.re),expj(z.im));
	  }
	  
	public static Complex expj(double theta)
	  { return new Complex(Math.cos(theta), Math.sin(theta));
	  }

	  public static Complex multiply(Complex z1, Complex z2){
		Complex result = new Complex (	z1.re*z2.re - z1.im*z2.im,
										z1.re*z2.im + z1.im*z2.re);
		return result;
	  }
	  
	  public static Complex changeSign(Complex z){
		return  new Complex(-z.re,-z.im);
	  }
	   
	  

	 public static Complex multiplyWithConstant(double a,Complex z){
		z.re *= a; z.im *= a;
		return z;
	}
	  
	  public static Complex divide(Complex z1, Complex z2){
		double mag = (z2.re * z2.re) + (z2.im * z2.im);
		Complex result = new Complex (((z1.re * z2.re) + (z1.im * z2.im)) / mag,
				((z1.im * z2.re) - (z1.re * z2.im)) / mag);
		return result;
	  }
	  
	 public static Complex divideConstantBy(double a,Complex z){
		z.re = a/z.re; z.im = a/z.im;
		return z;
	}
	  
	  public static Complex sum(Complex z1, Complex z2){
		z1.re += z2.re;
		z1.im += z2.im;
		return z1;
	  }
	  
	  public static Complex sumConstant(double a, Complex z){
		z.re += a;
		z.im += a;
		return z;
	  }
	  
	  public static Complex subtractFromConstant(double a, Complex z){
		z.re = a-z.re;
		z.im = a-z.im;
		return z;
	  }
	  
	  
	  public static Complex subtract(Complex z1, Complex z2){
		z1.re -= z2.re;
		z1.im -= z2.im;
		return z1;
	  }

}
