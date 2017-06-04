//  Copyright (c) 2013 Rocco Romeo
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


using System;
using System.Diagnostics;

namespace MathXK
{

    public static partial class Math2
    {

        /// <summary>
        /// Returns a ValueTuple: (Sin(x + piMultiple*π), Cos(x + piMultiple*π))
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="piMultiple">The PI multiple</param>
        internal static (double Sin, double Cos) SinCos(double x, double piMultiple)
        {
            if ((double.IsNaN(x) || double.IsInfinity(x))
            || (double.IsNaN(piMultiple) || double.IsInfinity(piMultiple))) {
                Policies.ReportDomainError("SinCos(x: {0}, piMultiple: {1}): Requires finite x, piMultiple", x, piMultiple);
                return (double.NaN, double.NaN);
            }

            // reduce multiple to (-2, 2)
            if (Math.Abs(piMultiple) >= 2)
                piMultiple = Math2.Mod(piMultiple, 2);

            double sinX = Math.Sin(x);
            double cosX = Math.Cos(x);

            if (piMultiple == 0)
                return (sinX, cosX);

            double sinY, cosY;

            double absMult = Math.Abs(piMultiple);
            int quadrant = (int)(absMult * 2);
            absMult -= 0.5 * quadrant; // absMult in [0, 0.5]

            double sinA, cosA;
            if (absMult > 0.25) {
                sinA = Math.Cos((0.5 - absMult) * Math.PI);
                cosA = Math.Sin((0.5 - absMult) * Math.PI);
            } else {
                sinA = Math.Sin(absMult * Math.PI);
                cosA = Math.Cos(absMult * Math.PI);
            }

            switch (quadrant) {
            case 0:
                sinY = sinA;
                cosY = cosA;
                break;
            case 1:
                // Sin(PI/2 + x) = Cos(x), Cos(PI/2 + x) = -Sin(x)
                sinY = cosA;
                cosY = -sinA;
                break;
            case 2:
                // Sin(PI + x) = -Sin(x), Cos(PI + x) = -Cos(x)
                sinY = -sinA;
                cosY = -cosA;
                break;
            case 3:
                // Sin(3*PI/2 + x) = -Cos(x), Cos(3*PI/2 + x) = Sin(x)
                sinY = -cosA;
                cosY = sinA;
                break;
            default:
                Policies.ReportEvaluationError("SinCosPI(piMultiple: {0}): Evaluation Error", piMultiple);
                return (double.NaN, double.NaN);
            }

            // sin(-x) == -sin(x), cos(-x) == cos(x)
            if (piMultiple < 0)
                sinY = -sinY;

            double sin = sinX * cosY + cosX * sinY;
            double cos = cosX * cosY - sinX * sinY;


            return (sin, cos);

        }

    }


} // namespaces

