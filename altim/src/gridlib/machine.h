// $Id: machine.h,v 1.1.1.1 1999/02/23 17:51:13 remko Exp $
//
// This is a machine dependent file, no promises that it will
// work in general, it does on a hp9000/735 under HP/UX. In case
// you ever go to another machine, make sure that:
//
//   int2   --> 2 byte signed integer
//   int4   --> 4 byte signed integer
//   real4  --> 4 byte floating point number
//   real8  --> 8 byte floating point number
//   real16 --> 16 byte floating point number
//
// as is done in the #defines just below. Some caveats:
//
//   - On a PC using Borland C/C++ 3.1 I noticed that an int2 is 
//     defined as an int. The best way to find out is to use the
//     debugger.
//
#define int2   short
#define int4   int
#define real4  float
#define real8  double
#define real16 long double
//
// Then for booleans
//
#define TRUE   1
#define FALSE  0
enum boolean {False=0, True=1}; 
//
// E. Schrama, schrama@geo.tudelft.nl   17/12/96   (C) TUD 1996
