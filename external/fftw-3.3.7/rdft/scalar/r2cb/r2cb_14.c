/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Sun Oct 29 08:18:36 EDT 2017 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_r2cb.native -fma -compact -variables 4 -pipeline-latency 4 -sign 1 -n 14 -name r2cb_14 -include rdft/scalar/r2cb.h */

/*
 * This function contains 62 FP additions, 44 FP multiplications,
 * (or, 18 additions, 0 multiplications, 44 fused multiply/add),
 * 46 stack variables, 7 constants, and 28 memory accesses
 */
#include "rdft/scalar/r2cb.h"

static void r2cb_14(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP1_949855824, +1.949855824363647214036263365987862434465571601);
     DK(KP1_801937735, +1.801937735804838252472204639014890102331838324);
     DK(KP692021471, +0.692021471630095869627814897002069140197260599);
     DK(KP801937735, +0.801937735804838252472204639014890102331838324);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     DK(KP356895867, +0.356895867892209443894399510021300583399127187);
     DK(KP554958132, +0.554958132087371191422194871006410481067288862);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ovs, R1 = R1 + ovs, Cr = Cr + ivs, Ci = Ci + ivs, MAKE_VOLATILE_STRIDE(56, rs), MAKE_VOLATILE_STRIDE(56, csr), MAKE_VOLATILE_STRIDE(56, csi)) {
	       E T3, Te, To, TK, Tu, TM, Tr, TL, Tv, TA, TX, TS, TN, TF, T6;
	       E Tf, Tc, Th, T9, Tg, Tj, Tx, TU, TP, TH, TC, T1, T2, Td, Ti;
	       T1 = Cr[0];
	       T2 = Cr[WS(csr, 7)];
	       T3 = T1 - T2;
	       Te = T1 + T2;
	       {
		    E Tm, Tn, T4, T5;
		    Tm = Ci[WS(csi, 4)];
		    Tn = Ci[WS(csi, 3)];
		    To = Tm - Tn;
		    TK = Tm + Tn;
		    {
			 E Ts, Tt, Tp, Tq;
			 Ts = Ci[WS(csi, 6)];
			 Tt = Ci[WS(csi, 1)];
			 Tu = Ts - Tt;
			 TM = Ts + Tt;
			 Tp = Ci[WS(csi, 2)];
			 Tq = Ci[WS(csi, 5)];
			 Tr = Tp - Tq;
			 TL = Tp + Tq;
		    }
		    Tv = FMA(KP554958132, Tu, Tr);
		    TA = FMA(KP554958132, To, Tu);
		    TX = FNMS(KP554958132, TL, TK);
		    TS = FMA(KP554958132, TK, TM);
		    TN = FMA(KP554958132, TM, TL);
		    TF = FNMS(KP554958132, Tr, To);
		    T4 = Cr[WS(csr, 2)];
		    T5 = Cr[WS(csr, 5)];
		    T6 = T4 - T5;
		    Tf = T4 + T5;
		    {
			 E Ta, Tb, T7, T8;
			 Ta = Cr[WS(csr, 6)];
			 Tb = Cr[WS(csr, 1)];
			 Tc = Ta - Tb;
			 Th = Ta + Tb;
			 T7 = Cr[WS(csr, 4)];
			 T8 = Cr[WS(csr, 3)];
			 T9 = T7 - T8;
			 Tg = T7 + T8;
		    }
		    Tj = FNMS(KP356895867, Tg, Tf);
		    Tx = FNMS(KP356895867, Tf, Th);
		    TU = FNMS(KP356895867, Tc, T9);
		    TP = FNMS(KP356895867, T6, Tc);
		    TH = FNMS(KP356895867, T9, T6);
		    TC = FNMS(KP356895867, Th, Tg);
	       }
	       Td = T6 + T9 + Tc;
	       R1[WS(rs, 3)] = FMA(KP2_000000000, Td, T3);
	       Ti = Tf + Tg + Th;
	       R0[0] = FMA(KP2_000000000, Ti, Te);
	       {
		    E Tw, Tl, Tk, TY, TW, TV;
		    Tw = FMA(KP801937735, Tv, To);
		    Tk = FNMS(KP692021471, Tj, Th);
		    Tl = FNMS(KP1_801937735, Tk, Te);
		    R0[WS(rs, 4)] = FNMS(KP1_949855824, Tw, Tl);
		    R0[WS(rs, 3)] = FMA(KP1_949855824, Tw, Tl);
		    TY = FNMS(KP801937735, TX, TM);
		    TV = FNMS(KP692021471, TU, T6);
		    TW = FNMS(KP1_801937735, TV, T3);
		    R1[WS(rs, 1)] = FNMS(KP1_949855824, TY, TW);
		    R1[WS(rs, 5)] = FMA(KP1_949855824, TY, TW);
	       }
	       {
		    E TB, Tz, Ty, TO, TJ, TI;
		    TB = FNMS(KP801937735, TA, Tr);
		    Ty = FNMS(KP692021471, Tx, Tg);
		    Tz = FNMS(KP1_801937735, Ty, Te);
		    R0[WS(rs, 1)] = FNMS(KP1_949855824, TB, Tz);
		    R0[WS(rs, 6)] = FMA(KP1_949855824, TB, Tz);
		    TO = FMA(KP801937735, TN, TK);
		    TI = FNMS(KP692021471, TH, Tc);
		    TJ = FNMS(KP1_801937735, TI, T3);
		    R1[0] = FNMS(KP1_949855824, TO, TJ);
		    R1[WS(rs, 6)] = FMA(KP1_949855824, TO, TJ);
	       }
	       {
		    E TT, TR, TQ, TG, TE, TD;
		    TT = FNMS(KP801937735, TS, TL);
		    TQ = FNMS(KP692021471, TP, T9);
		    TR = FNMS(KP1_801937735, TQ, T3);
		    R1[WS(rs, 4)] = FNMS(KP1_949855824, TT, TR);
		    R1[WS(rs, 2)] = FMA(KP1_949855824, TT, TR);
		    TG = FNMS(KP801937735, TF, Tu);
		    TD = FNMS(KP692021471, TC, Tf);
		    TE = FNMS(KP1_801937735, TD, Te);
		    R0[WS(rs, 5)] = FNMS(KP1_949855824, TG, TE);
		    R0[WS(rs, 2)] = FMA(KP1_949855824, TG, TE);
	       }
	  }
     }
}

static const kr2c_desc desc = { 14, "r2cb_14", {18, 0, 44, 0}, &GENUS };

void X(codelet_r2cb_14) (planner *p) {
     X(kr2c_register) (p, r2cb_14, &desc);
}

#else

/* Generated by: ../../../genfft/gen_r2cb.native -compact -variables 4 -pipeline-latency 4 -sign 1 -n 14 -name r2cb_14 -include rdft/scalar/r2cb.h */

/*
 * This function contains 62 FP additions, 38 FP multiplications,
 * (or, 36 additions, 12 multiplications, 26 fused multiply/add),
 * 28 stack variables, 7 constants, and 28 memory accesses
 */
#include "rdft/scalar/r2cb.h"

static void r2cb_14(R *R0, R *R1, R *Cr, R *Ci, stride rs, stride csr, stride csi, INT v, INT ivs, INT ovs)
{
     DK(KP1_801937735, +1.801937735804838252472204639014890102331838324);
     DK(KP445041867, +0.445041867912628808577805128993589518932711138);
     DK(KP1_246979603, +1.246979603717467061050009768008479621264549462);
     DK(KP867767478, +0.867767478235116240951536665696717509219981456);
     DK(KP1_949855824, +1.949855824363647214036263365987862434465571601);
     DK(KP1_563662964, +1.563662964936059617416889053348115500464669037);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, R0 = R0 + ovs, R1 = R1 + ovs, Cr = Cr + ivs, Ci = Ci + ivs, MAKE_VOLATILE_STRIDE(56, rs), MAKE_VOLATILE_STRIDE(56, csr), MAKE_VOLATILE_STRIDE(56, csi)) {
	       E T3, Td, T6, Te, Tq, Tz, Tn, Ty, Tc, Tg, Tk, Tx, T9, Tf, T1;
	       E T2;
	       T1 = Cr[0];
	       T2 = Cr[WS(csr, 7)];
	       T3 = T1 - T2;
	       Td = T1 + T2;
	       {
		    E T4, T5, To, Tp;
		    T4 = Cr[WS(csr, 2)];
		    T5 = Cr[WS(csr, 5)];
		    T6 = T4 - T5;
		    Te = T4 + T5;
		    To = Ci[WS(csi, 2)];
		    Tp = Ci[WS(csi, 5)];
		    Tq = To - Tp;
		    Tz = To + Tp;
	       }
	       {
		    E Tl, Tm, Ta, Tb;
		    Tl = Ci[WS(csi, 6)];
		    Tm = Ci[WS(csi, 1)];
		    Tn = Tl - Tm;
		    Ty = Tl + Tm;
		    Ta = Cr[WS(csr, 6)];
		    Tb = Cr[WS(csr, 1)];
		    Tc = Ta - Tb;
		    Tg = Ta + Tb;
	       }
	       {
		    E Ti, Tj, T7, T8;
		    Ti = Ci[WS(csi, 4)];
		    Tj = Ci[WS(csi, 3)];
		    Tk = Ti - Tj;
		    Tx = Ti + Tj;
		    T7 = Cr[WS(csr, 4)];
		    T8 = Cr[WS(csr, 3)];
		    T9 = T7 - T8;
		    Tf = T7 + T8;
	       }
	       R1[WS(rs, 3)] = FMA(KP2_000000000, T6 + T9 + Tc, T3);
	       R0[0] = FMA(KP2_000000000, Te + Tf + Tg, Td);
	       {
		    E Tr, Th, TE, TD;
		    Tr = FNMS(KP1_949855824, Tn, KP1_563662964 * Tk) - (KP867767478 * Tq);
		    Th = FMA(KP1_246979603, Tf, Td) + FNMA(KP445041867, Tg, KP1_801937735 * Te);
		    R0[WS(rs, 2)] = Th - Tr;
		    R0[WS(rs, 5)] = Th + Tr;
		    TE = FMA(KP867767478, Tx, KP1_563662964 * Ty) - (KP1_949855824 * Tz);
		    TD = FMA(KP1_246979603, Tc, T3) + FNMA(KP1_801937735, T9, KP445041867 * T6);
		    R1[WS(rs, 2)] = TD - TE;
		    R1[WS(rs, 4)] = TD + TE;
	       }
	       {
		    E Tt, Ts, TA, Tw;
		    Tt = FMA(KP867767478, Tk, KP1_563662964 * Tn) - (KP1_949855824 * Tq);
		    Ts = FMA(KP1_246979603, Tg, Td) + FNMA(KP1_801937735, Tf, KP445041867 * Te);
		    R0[WS(rs, 6)] = Ts - Tt;
		    R0[WS(rs, 1)] = Ts + Tt;
		    TA = FNMS(KP1_949855824, Ty, KP1_563662964 * Tx) - (KP867767478 * Tz);
		    Tw = FMA(KP1_246979603, T9, T3) + FNMA(KP445041867, Tc, KP1_801937735 * T6);
		    R1[WS(rs, 5)] = Tw - TA;
		    R1[WS(rs, 1)] = Tw + TA;
	       }
	       {
		    E TC, TB, Tv, Tu;
		    TC = FMA(KP1_563662964, Tz, KP1_949855824 * Tx) + (KP867767478 * Ty);
		    TB = FMA(KP1_246979603, T6, T3) + FNMA(KP1_801937735, Tc, KP445041867 * T9);
		    R1[0] = TB - TC;
		    R1[WS(rs, 6)] = TB + TC;
		    Tv = FMA(KP1_563662964, Tq, KP1_949855824 * Tk) + (KP867767478 * Tn);
		    Tu = FMA(KP1_246979603, Te, Td) + FNMA(KP1_801937735, Tg, KP445041867 * Tf);
		    R0[WS(rs, 4)] = Tu - Tv;
		    R0[WS(rs, 3)] = Tu + Tv;
	       }
	  }
     }
}

static const kr2c_desc desc = { 14, "r2cb_14", {36, 12, 26, 0}, &GENUS };

void X(codelet_r2cb_14) (planner *p) {
     X(kr2c_register) (p, r2cb_14, &desc);
}

#endif
