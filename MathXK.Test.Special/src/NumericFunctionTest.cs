//  (C) Copyright Rocco Romeo 2013.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

using System;
using System.Collections.Generic;
using MathXK.Test.Functions;
using Microsoft.VisualStudio.TestTools.UnitTesting;


namespace MathXK.Test
{
    public class TestCase {

        public static TestCase<T1> Create<T1>(T1 param1,double result)
        {
            return new TestCase<T1>(param1, result);
        }

        public static TestCase<T1, T2> Create<T1, T2>(T1 param1, T2 param2, double result)
        {
            return new TestCase<T1, T2>(param1, param2, result);
        }

        public static TestCase<T1, T2, T3> Create<T1, T2, T3>(T1 param1, T2 param2, T3 param3, double result)
        {
            return new TestCase<T1, T2, T3>(param1, param2, param3, result);
        }

        public static TestCase<T1, T2, T3, T4> Create<T1, T2, T3, T4>(T1 param1, T2 param2, T3 param3, T4 param4, double result)
        {
            return new TestCase<T1, T2, T3, T4>(param1, param2, param3, param4, result);
        }

    }

    public class TestCase<T1>
    {
        public T1 Param1;
        public double Result;

        public TestCase(T1 param1, double result)
        {
            Param1 = param1;
            Result = result;
        }

        public override string ToString()
        {
            return string.Format("({0}) = {1:E16}", Param1, Result);
        }


    }

    public class TestCase<T1, T2>
    {
        public T1 Param1;
        public T2 Param2;
        public double Result;

        public TestCase(T1 param1, T2 param2, double result)
        {
            Param1 = param1;
            Param2 = param2;
            Result = result;
        }

        public override string ToString()
        {
            return string.Format("({0}, {1}) = {2:E16}", Param1, Param2, Result);
        }
    }

    public class TestCase<T1, T2, T3>
    {
        public T1 Param1;
        public T2 Param2;
        public T3 Param3;
        public double Result;

        public TestCase(T1 param1, T2 param2, T3 param3, double result)
        {
            Param1 = param1;
            Param2 = param2;
            Param3 = param3;
            Result = result;
        }

        public override string ToString()
        {
            return string.Format("({0}, {1}, {2}) = {3:E16}", Param1, Param2, Param3, Result);
        }
    }

    public class TestCase<T1, T2, T3, T4>
    {
        public T1 Param1;
        public T2 Param2;
        public T3 Param3;
        public T4 Param4;
        public double Result;

        public TestCase(T1 param1, T2 param2, T3 param3, T4 param4, double result)
        {
            Param1 = param1;
            Param2 = param2;
            Param3 = param3;
            Param4 = param4;
            Result = result;
        }

        public override string ToString()
        {
            return string.Format("({0}, {1}, {2}, {3}) = {4:E16}", Param1, Param2, Param3, Param4, Result);
        }
    }


    //=======================
    // TestCaseSet
    //========================

    public class TestCaseSet<T1>
    {
        public IEnumerable<TestCase<T1>> DataSet;
        public string DataSetName;
        public int Ulps;
        public int AveUlps;

        public TestCaseSet(IEnumerable<TestCase<T1>> data, string name, int ulps, int aveUlps)
        {
            DataSet = data;
            DataSetName = name;
            Ulps = ulps;
            AveUlps = aveUlps;

        }
        public TestCaseSet(IEnumerable<TestCase<T1>> data, string name, int ulps) : this(data, name, ulps, ulps) { }
    }

    public class TestCaseSet<T1, T2>
    {
        public IEnumerable<TestCase<T1, T2>> DataSet;
        public string DataSetName;
        public int Ulps;
        public int AveUlps;

        public TestCaseSet(IEnumerable<TestCase<T1, T2>> data, string name, int ulps, int aveUlps)
        {
            DataSet = data;
            DataSetName = name;
            Ulps = ulps;
            AveUlps = aveUlps;
        }
        public TestCaseSet(IEnumerable<TestCase<T1, T2>> data, string name, int ulps) : this(data, name, ulps, ulps) { }
    }

    public class TestCaseSet<T1, T2, T3>
    {
        public IEnumerable<TestCase<T1, T2, T3>> DataSet;
        public string DataSetName;
        public int Ulps;
        public int AveUlps;

        public TestCaseSet(IEnumerable<TestCase<T1, T2, T3>> data, string name, int ulps, int aveUlps)
        {
            DataSet = data;
            DataSetName = name;
            Ulps = ulps;
            AveUlps = aveUlps;
        }

        public TestCaseSet(IEnumerable<TestCase<T1, T2, T3>> data, string name, int ulps) : this(data, name, ulps, ulps) { }
    }

    public class TestCaseSet<T1, T2, T3, T4>
    {
        public IEnumerable<TestCase<T1, T2, T3, T4>> DataSet;
        public string DataSetName;
        public int Ulps;
        public int AveUlps;

        public TestCaseSet(IEnumerable<TestCase<T1, T2, T3, T4>> data, string name, int ulps, int aveUlps)
        {
            DataSet = data;
            DataSetName = name;
            Ulps = ulps;
            AveUlps = aveUlps;
        }

        public TestCaseSet(IEnumerable<TestCase<T1, T2, T3, T4>> data, string name, int ulps) : this(data, name, ulps, ulps) { }
    }


    class NumericFunctionTest
    {
        private static bool CatchExceptions()
        {
#if THROW
            return true;
#else
            return false;
#endif
        }

        public static double Distance(double a, double b)
        {
            if (double.IsNaN(a) || double.IsNaN(b))
                return (double.IsNaN(a) && double.IsNaN(b)) ? 0 : double.MaxValue;

            if (double.IsInfinity(a) || double.IsInfinity(b))
                return (double.IsInfinity(a) && double.IsInfinity(b) && Math.Sign(a) == Math.Sign(b)) ? 0 : double.MaxValue;

            double result = Math2.FloatDistance(a, b);
            if (double.IsNaN(result))
                return double.MaxValue;

            return result;
        }


        public static void RunSet<T1>(Func<T1, double> f, string funcName, IEnumerable<TestCaseSet<T1>> testCases)
        {

            foreach(var testCase in testCases) {

                var dataSet = testCase.DataSet;
                int maxUlps = testCase.Ulps;
                string desc = testCase.DataSetName;
                int maxAveUlps = testCase.AveUlps;

                double ulpsSum = 0;

                int n = 0;
                foreach (var item in dataSet) {
                    T1 param1 = item.Param1;
                    double expected = item.Result;
                    double calculated;

                    try {
                        calculated = f(param1);
                    } catch (Exception e) {
                        if (!CatchExceptions())
                            throw;

                        if ((e is Exceptions.DomainException)
                            || (e is Exceptions.PoleException)
                            || (e is Exceptions.ConvergenceException)
                            || (e is Exceptions.NotImplementedException)
                            || (e is Exceptions.RootNotFoundException)) {
                            calculated = double.NaN;
                        } else {
                            throw;
                        }
                    }

                    double ulps = Distance(expected, calculated);
                    if (Math.Abs(ulps) > maxUlps) {
                        string output = $"{funcName}{item}; got {calculated:E16}; MaxUlps = {maxUlps}; ulps = {ulps}; DataSet = {desc}";
                        throw new AssertFailedException(output);
                    }

                    ulpsSum += Math.Abs(ulps);
                    n++;
                }

                double aveUlps = (n == 0) ? 0 : ulpsSum / n;
                if (aveUlps > maxAveUlps) {
                    string output = $"{funcName} Ave Ulps = {aveUlps}; MaxUlps = {maxAveUlps}; DataSet = {desc}";
                    throw new AssertFailedException(output);
                }

            }

        }

        public static void RunSet<T1, T2>(Func<T1, T2, double> f, string funcName, IEnumerable<TestCaseSet<T1, T2>> testCases)
        {

            foreach (var testCase in testCases) {

                var dataSet = testCase.DataSet;
                int maxUlps = testCase.Ulps;
                string desc = testCase.DataSetName;
                int maxAveUlps = testCase.AveUlps;

                double ulpsSum = 0;

                int n = 0;
                foreach (var item in dataSet) {
                    T1 param1 = item.Param1;
                    T2 param2 = item.Param2;
                    double expected = item.Result;

                    double calculated;

                    try {
                        calculated = f(param1, param2);
                    } catch (Exception e) {
                        if (!CatchExceptions())
                            throw;

                        if ((e is Exceptions.DomainException)
                            || (e is Exceptions.PoleException)
                            || (e is Exceptions.ConvergenceException)
                            || (e is Exceptions.NotImplementedException)
                            || (e is Exceptions.RootNotFoundException)) {
                            calculated = double.NaN;
                        } else {
                            throw;
                        }
                    }

                    double ulps = Distance(expected, calculated);
                    if (Math.Abs(ulps) > maxUlps) {
                        string output = $"{funcName}{item}; got {calculated:E16}; MaxUlps = {maxUlps}; ulps = {ulps}; DataSet = {desc}";
                        throw new AssertFailedException(output);
                    }

                    ulpsSum += Math.Abs(ulps);
                    n++;
                }

                double aveUlps = (n == 0) ? 0 : ulpsSum / n;
                if (aveUlps > maxAveUlps) {
                    string output = $"{funcName} Ave Ulps = {aveUlps}; MaxUlps = {maxAveUlps}; DataSet = {desc}";
                    throw new AssertFailedException(output);
                }

            }

        }

        public static void RunSet<T1, T2, T3>(Func<T1, T2, T3, double> f, string funcName, IEnumerable<TestCaseSet<T1, T2, T3>> testCases)
        {

            foreach (var testCase in testCases) {

                var dataSet = testCase.DataSet;
                int maxUlps = testCase.Ulps;
                string desc = testCase.DataSetName;
                int maxAveUlps = testCase.AveUlps;

                double ulpsSum = 0;

                int n = 0;
                foreach (var item in dataSet) {
                    T1 param1 = item.Param1;
                    T2 param2 = item.Param2;
                    T3 param3 = item.Param3;
                    double expected = item.Result;

                    double calculated;

                    try {
                        calculated = f(param1, param2, param3);
                    } catch (Exception e) {
                        if (!CatchExceptions())
                            throw;

                        if ((e is Exceptions.DomainException)
                            || (e is Exceptions.PoleException)
                            || (e is Exceptions.ConvergenceException)
                            || (e is Exceptions.NotImplementedException)
                            || (e is Exceptions.RootNotFoundException)) {
                            calculated = double.NaN;
                        } else {
                            throw;
                        }
                    }


                    double ulps = Distance(expected, calculated);
                    if (Math.Abs(ulps) > maxUlps) {
                        string output = $"{funcName}{item}; got {calculated:E16}; MaxUlps = {maxUlps}; ulps = {ulps}; DataSet = {desc}";
                        throw new AssertFailedException(output);
                    }
                    ulpsSum += Math.Abs(ulps);
                    n++;
                }

                double aveUlps = (n == 0) ? 0 : ulpsSum / n;
                if (aveUlps > maxAveUlps) {
                    string output = $"{funcName} Ave Ulps = {aveUlps}; MaxUlps = {maxAveUlps}; DataSet = {desc}";
                    throw new AssertFailedException(output);
                }

            }

        }

        public static void RunSet<T1, T2, T3, T4>(Func<T1, T2, T3, T4, double> f, string funcName, IEnumerable<TestCaseSet<T1, T2, T3, T4>> testCases)
        {

            foreach (var testCase in testCases) {

                var dataSet = testCase.DataSet;
                int maxUlps = testCase.Ulps;
                string desc = testCase.DataSetName;
                int maxAveUlps = testCase.AveUlps;

                double ulpsSum = 0;

                int n = 0;
                foreach (var item in dataSet) {
                    T1 param1 = item.Param1;
                    T2 param2 = item.Param2;
                    T3 param3 = item.Param3;
                    T4 param4 = item.Param4;
                    double expected = item.Result;

                    double calculated;

                    try {
                        calculated = f(param1, param2, param3, param4);
                    } catch (Exception e) {
                        if (!CatchExceptions())
                            throw;

                        if ((e is Exceptions.DomainException)
                            || (e is Exceptions.PoleException)
                            || (e is Exceptions.ConvergenceException)
                            || (e is Exceptions.NotImplementedException)
                            || (e is Exceptions.RootNotFoundException)) {
                            calculated = double.NaN;
                        } else {
                            throw;
                        }
                    }

                    double ulps = Distance(expected, calculated);
                    if (Math.Abs(ulps) > maxUlps) {
                        string output = $"{funcName}{item}; got {calculated:E16}; MaxUlps = {maxUlps}; ulps = {ulps}; DataSet = {desc}";
                        throw new AssertFailedException(output);
                    }
                    ulpsSum += Math.Abs(ulps);
                    n++;
                }

                double aveUlps = (n == 0) ? 0 : ulpsSum / n;
                if (aveUlps > maxAveUlps) {
                    string output = $"{funcName} Ave Ulps = {aveUlps}; MaxUlps = {maxAveUlps}; DataSet = {desc}";
                    throw new AssertFailedException(output);
                }

            }

        }

    }



}
