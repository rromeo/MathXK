//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Based on the work of:
//      (C) Copyright John Maddock & Paul A. Bristow 2007, Boost Software License v1.0


// ===========================================================
// The default error handling behavior is as follows:
// DEBUG - Report a dignostic message to debug and return NaN
// RELEASE - return NaN
//
// if you prefer to throw an exception instead
// define THROW in the Visual Studio Project Properties
// ============================================================

#if THROW
#define MATHXK_THROW_DOMAIN_ERROR
#define MATHXK_THROW_POLE_ERROR
#define MATHXK_THROW_NOTIMPLEMENTED_ERROR
#define MATHXK_THROW_EVALUATION_ERROR
#define MATHXK_THROW_CONVERGENCE_ERROR
#define MATHXK_THROW_ROOTNOTFOUND_ERROR

#endif

#if DEBUG
#define MATHXK_REPORT_DOMAIN_ERROR
#define MATHXK_REPORT_POLE_ERROR
#define MATHXK_REPORT_NOTIMPLEMENTED_ERROR
#define MATHXK_REPORT_EVALUATION_ERROR
#define MATHXK_REPORT_CONVERGENCE_ERROR
#define MATHXK_REPORT_ROOTNOTFOUND_ERROR

#endif




using System;
using System.Diagnostics;

namespace MathXK
{

    internal static class Policies
    {

        /// <summary>
        /// Maximum number of iterations permitted in series sum calculations before raising 
        /// a Convergence Error
        /// <para>Currently = 1000</para>
        /// </summary>
        public const int MaxSeriesIterations = 1000;

        /// <summary>
        /// Maximum number of iterations permitted in root finding algorithms before raising 
        /// a Root Not Found Error
        /// </summary>
        public const int MaxRootIterations = 200;

        /// <summary>
        /// Maximum number of iterations permitted in continued fractions algorithms before raising 
        /// a Convergence Exception
        /// </summary>
        public const int MaxFractionIterations = 200;

        /// <summary>
        /// Series convergence tolerance, unless otherwise specified
        /// </summary>
        public const double SeriesTolerance = 1 * DoubleLimits.MachineEpsilon;

        // ===============================================================================
        // Error/Warning Reporting
        // ===============================================================================

        /// <summary>
        /// Reported when more or more arguments are outside the defined range of the function
        /// </summary>
        /// <param name="message">The reason for the report as a format string</param>
        /// <param name="args">Additional arguments</param>
        [Conditional("DEBUG"), Conditional("THROW")]
        public static void ReportDomainError(String message, params object[] args)
        {
#if MATHXK_REPORT_DOMAIN_ERROR
            Debug.WriteLine(message, args);
#endif

#if MATHXK_THROW_DOMAIN_ERROR
            throw new Exceptions.DomainException(message, args);
#endif
        }

        /// <summary>
        /// Reported when more or more arguments would cause the function to be evaluated at a pole
        /// </summary>
        /// <param name="message">The reason for the report as a format string</param>
        /// <param name="args">Additional arguments</param>
        [Conditional("DEBUG"), Conditional("THROW")]
        public static void ReportPoleError(String message, params object[] args)
        {
#if MATHXK_REPORT_POLE_ERROR
            Debug.WriteLine(message, args);
#endif

#if MATHXK_THROW_POLE_ERROR
            throw new Exceptions.PoleException(message, args);
#endif
        }



        /// <summary>
        /// Reported when a function has not implemented entirely or for a specific parameter range
        /// </summary>
        /// <param name="message">The reason for the report as a format string</param>
        /// <param name="args">Additional arguments</param>
        [Conditional("DEBUG"), Conditional("THROW")]
        public static void ReportNotImplementedError(String message, params object[] args)
        {
#if MATHXK_REPORT_NOTIMPLEMENTED_ERROR
            Debug.WriteLine(message, args);
#endif

#if MATHXK_THROW_NOTIMPLEMENTED_ERROR
            throw new Exceptions.NotImplementedException(message, args);
#endif
        }


        /// <summary>
        /// Reported when the result of the function should have been computed properly, but was not
        /// </summary>
        /// <param name="message">The reason for the report as a format string</param>
        /// <param name="args">Additional arguments</param>
        [Conditional("DEBUG"), Conditional("THROW")]
        public static void ReportEvaluationError(String message, params object[] args)
        {
#if MATHXK_REPORT_EVALUATION_ERROR
            Debug.WriteLine(message, args);
#endif

#if MATHXK_THROW_EVALUATION_ERROR
            throw new Exceptions.EvaluationException(message, args);
#endif
        }

        /// <summary>
        /// Reported when a series has not converged within the specified number of iterations
        /// </summary>
        /// <param name="message">The reason for the report as a format string</param>
        /// <param name="args">Additional arguments</param>
        [Conditional("DEBUG"), Conditional("THROW")]
        public static void ReportConvergenceError(String message, params object[] args)
        {
#if MATHXK_REPORT_CONVERGENCE_ERROR
            Debug.WriteLine(message, args);
#endif

#if MATHXK_THROW_CONVERGENCE_ERROR
            throw new Exceptions.ConvergenceException(message, args);
#endif
        }

        /// <summary>
        /// Reported when a root(zero) has not been found within the specified number of iterations or tolerance
        /// </summary>
        /// <param name="message">The reason for the report as a format string</param>
        /// <param name="args">Additional arguments</param>
        [Conditional("DEBUG"), Conditional("THROW")]
        public static void ReportRootNotFoundError(String message, params object[] args)
        {
#if MATHXK_REPORT_ROOTNOTFOUND_ERROR
            Debug.WriteLine(message, args);
#endif

#if MATHXK_THROW_ROOTNOTFOUND_ERROR
            throw new Exceptions.RootNotFoundException(message, args);
#endif
        }



        // ===============================================================================
        // Error Handling
        // ===============================================================================


        [Conditional("DEBUG"), Conditional("EXTRA_DEBUG")]
        internal static void WriteLine(string message, params object[] args)
        {
            Debug.WriteLine(message, args);
        }

    }


} // namespaces


