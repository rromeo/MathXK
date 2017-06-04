//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Based on the work of:
//      Copyright (c) John Maddock & Paul A. Bristow 2007, Boost Software License v1.0


namespace MathXK.Exceptions
{

    ///<summary>Raised when more or more arguments are outside the defined range of the function</summary>
    public class DomainException: System.ArgumentException
    {
        /// <summary>
        /// Creates a Domain Exception with the given message
        /// </summary>
        /// <param name="message">The reason for the exception</param>
        public DomainException(string message) : base(message) { }

        /// <summary>
        ///  Creates a Domain Exception with the given message and arguments
        /// </summary>
        /// <param name="message">The reason for the exception as a format string</param>
        /// <param name="args">Additional arguments</param>
        public DomainException(string message, params object[] args) : base(string.Format(message, args)) { }


        /// <summary>
        ///  Creates a Domain Exception with the given message and inner exception
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="inner">The inner exception</param>
        public DomainException(string message, System.Exception inner) : base(message, inner) { }
    }

    ///<summary>Raised when more or more arguments would cause the function to be evaluated at a pole</summary>
    public class PoleException: System.ArgumentException
    {
        /// <summary>
        /// Creates a Pole Exception with the given message
        /// </summary>
        /// <param name="message">The reason for the exception</param>
        public PoleException(string message) : base(message) { }

        /// <summary>
        ///  Creates a Pole Exception with the given message and arguments
        /// </summary>
        /// <param name="message">The reason for the exception as a format string</param>
        /// <param name="args">Additional arguments</param>
        public PoleException(string message, params object[] args) : base(string.Format(message, args)) { }

        /// <summary>
        ///  Creates a Pole Exception with the given message and inner exception
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="inner">The inner exception</param>
        public PoleException(string message, System.Exception inner) : base(message, inner) { }

    }


    ///<summary>Raised when the result of the function is outside the representable range of the floating point type used</summary>
    public class OverflowException: System.OverflowException
    {
        /// <summary>
        /// Creates an Overflow Exception with the given message
        /// </summary>
        /// <param name="message">The reason for the exception</param>
        public OverflowException(string message) : base(message) { }


        /// <summary>
        ///  Creates an Overflow Exception with the given message and arguments
        /// </summary>
        /// <param name="message">The reason for the exception as a format string</param>
        /// <param name="args">Additional arguments</param>
        public OverflowException(string message, params object[] args) : base(string.Format(message, args)) { }


        /// <summary>
        ///  Creates an Overflow Exception with the given message and inner exception
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="inner">The inner exception</param>
        public OverflowException(string message, System.Exception inner) : base(message, inner) { }
    }

    ///<summary>Raised when the result of the function is too small to be represented in the floating point type used</summary>
    public class UnderflowException: System.ArithmeticException
    {
        /// <summary>
        /// Creates an Underflow Exception with the given message
        /// </summary>
        /// <param name="message">The reason for the exception</param>
        public UnderflowException(string message) : base(message) { }

        /// <summary>
        ///  Creates an Undeflow Exception with the given message and arguments
        /// </summary>
        /// <param name="message">The reason for the exception as a format string</param>
        /// <param name="args">Additional arguments</param>
        public UnderflowException(string message, params object[] args) : base(string.Format(message, args)) { }


        /// <summary>
        ///  Creates an Undeflow Exception with the given message and inner exception
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="inner">The inner exception</param>
        public UnderflowException(string message, System.Exception inner) : base(message, inner) { }
    }

    ///<summary>Raised when the result of the function should have been computed properly, but was not</summary>
    public class EvaluationException: System.ArithmeticException
    {
        /// <summary>
        /// Creates an Evaluation Exception with the given message
        /// </summary>
        /// <param name="message">The reason for the exception</param>
        public EvaluationException(string message) : base(message) { }

        /// <summary>
        ///  Creates an Evaluation Exception with the given message and arguments
        /// </summary>
        /// <param name="message">The reason for the exception as a format string</param>
        /// <param name="args">Additional arguments</param>
        public EvaluationException(string message, params object[] args) : base(string.Format(message, args)) { }

        /// <summary>
        ///  Creates an Evaluation Exception with the given message and inner exception
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="inner">The inner exception</param>
        public EvaluationException(string message, System.Exception inner) : base(message, inner) { }

    }

    ///<summary>Raised when a function has not implemented entirely or for a specific parameter range</summary>
    public class NotImplementedException: System.NotImplementedException
    {
        /// <summary>
        /// Creates a Not Implemented Exception
        /// </summary>
        public NotImplementedException() : base() { }

        /// <summary>
        /// Creates a Not Implemented Exception with the given message
        /// </summary>
        /// <param name="message">The reason for the exception</param>
        public NotImplementedException(string message) : base(message) { }

        /// <summary>
        /// Creates a Not Implemented Exception with the given message and arguments
        /// </summary>
        /// <param name="message">The reason for the exception as a format string</param>
        /// <param name="args">Additional arguments</param>
        public NotImplementedException(string message, params object[] args) : base(string.Format(message, args)) { }


        /// <summary>
        ///  Creates a Not Implemented Exception with the given message and inner exception
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="inner">The inner exception</param>
        public NotImplementedException(string message, System.Exception inner) : base(message, inner) { }

    }

    ///<summary>Raised when a series has not converged within the specified number of iterations</summary>
    public class ConvergenceException: System.ArithmeticException
    {
        /// <summary>
        /// Creates a Convergence Exception with the given message
        /// </summary>
        /// <param name="message">The reason for the exception</param>
        public ConvergenceException(string message) : base(message) { }

        /// <summary>
        /// Creates a Convergence Exception with the given message and arguments
        /// </summary>
        /// <param name="message">The reason for the exception as a format string</param>
        /// <param name="args">Additional arguments</param>
        public ConvergenceException(string message, params object[] args) : base(string.Format(message, args)) { }


        /// <summary>
        ///  Creates a Convergence Exception with the given message and inner exception
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="inner">The inner exception</param>
        public ConvergenceException(string message, System.Exception inner) : base(message, inner) { }
    }

    ///<summary>
    /// Raised when a root(zero) has not been found within the specified number of iterations or tolerance
    /// </summary>
    public class RootNotFoundException: System.ArithmeticException
    {
        /// <summary>
        /// Creates a Root Not Found Exception with the given message
        /// </summary>
        /// <param name="message">The reason for the exception</param>
        public RootNotFoundException(string message) : base(message) { }

        /// <summary>
        /// Creates a Root Not Found Exception with the given message and arguments
        /// </summary>
        /// <param name="message">The reason for the exception as a format string</param>
        /// <param name="args">Additional arguments</param>
        public RootNotFoundException(string message, params object[] args) : base(string.Format(message, args)) { }

        /// <summary>
        ///  Creates a Root Not Found Exception with the given message and inner exception
        /// </summary>
        /// <param name="message">The message</param>
        /// <param name="inner">The inner exception</param>
        public RootNotFoundException(string message, System.Exception inner) : base(message, inner) { }

    }

}
