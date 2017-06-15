# Introduction
Math Expansion Kit (MathXK) is a free open-source C# special functions math library based on the C++ Boost Math Toolkit. 

Benefits:
* Implements over 90 high-quality functions including: factorials, gammas, betas, error functions, C floating point functions, exponential integrals, bessels, elliptical integrals, polynomials, root-finding routines, and more.
* Passes the standard Boost Math tests and is internally tested against an arbitrary precision library. 
* .NET Standard 1.1 library
* Licensed under the Boost Software License v1.0: usable for commercial or non-commercial purposes without having to release your source code.

# Project Status

The routines have been tested with Windows/Intel x64; testing is required on other platforms. 

The statistical distributions remain to be done.

# Requirements
**Using the library (binaries) requires:**
* .NET Core 1.0+

.NET Framework will be supported when .NET standard 2.0 is released, apparently in [Q3 2017](https://github.com/dotnet/core/blob/master/roadmap.md)

**Compiling the source code requires:**
* Visual Studio 2017 
* C#7+
* MS-Test V2

The source code uses the new .NET Core csproj format (which is unsupported by earlier versions of Visual Studio) and uses C#7 value tuples. 

# Getting Started

## A. Using the Binaries with NuGet:

The easiest way to get started is to get the `MathXK` package from NuGet.

**In Visual Studio:** 

Use "Manage NuGet Packages..." to search for `MathXK`. (Make sure to check the "include prerelease" checkbox). Install the latest version.

**In Visual Studio Package Manager Console:** 
```
Install-Package MathXK -version "0.5.0-alpha"
```

**On the command line:**
```
 cd your_project_folder 
 dotnet add package MathXK -v "0.5.0-alpha"
```

## B. Downloading the source code:
If you prefer to work directly with the project, you can clone it to your local machine:
```
git clone https://github.com/rromeo/MathXK
```

# Quck Start Documentation
The main namespaces in MathXK are:
1. **MathXK** - Contains the static class **Math2**, where all the special functions are located.
2. **MathXK.Roots** - Contains all the root finding routines.
3. **MathXK.Exceptions** - Contains all the exception classes. 

A simple C# program that demonstrates the use of Math Extensions is shown below:
```C#
using System;

using MathXK; // Main namespace. Special functions are in Math2
using MathXK.Roots; // Root-finding namespace. Routines are in RootFinder

namespace MathXKExample
{
    class Program
    {
        static void Main(string[] args)
        {
            // Compute Gamma(0.5)
            var result = Math2.Tgamma(0.5);
            Console.WriteLine($"Gamma(0.5) = {result}");

            // Compute Beta(2, 1)
            result = Math2.Beta(2, 1);
            Console.WriteLine($"Beta(2, 1) = {result}");
            
            // Compute Asin(0.5)
            result = SimpleArcSin(0.5);
            Console.WriteLine($"ArcSin(0.5) = {result}");

            // Wait for input to close the console
            Console.Write("Press any key...");
            Console.ReadKey();
        }
        
        /// <summary>Arc Sin using root-finding</summary>
        static double SimpleArcSin(double radians) 
        {
            // the function and its derivative
            Func<double, (double, double)> f = (double x) => {
                var function = Math.Sin(x) - radians; 
                var derivative = Math.Cos(x);
                return (function, derivative);
            };

            var root = RootFinder.NewtonSafe(f, guess: 0, min: -1, max: 1);
            if ( root == null )
                return double.NaN; // Bad calling parameters
            if ( !root.Success )
                return double.NaN; // Root not found
            return root.SolutionX; // Return the root
        }
    }
}

```

### Error Handling
The default behavior of the library is to return `NaN` for function errors (similar to most System.Math functions). Only constructors and property setters throw exceptions.
 
If you wish the library to throw exceptions instead of returning `NaN`, define `THROW` in the Conditional Compilation Settings for the MathXK project and the unit test project.

# License
The Lift Math source code is licensed under the Boost Software License v1.0. 

"The Boost license permits the creation of derivative works for commercial or non-commercial use with no legal requirement to release your source code."  [Boost](http://www.boost.org/users/license.html)

For more information, see:
* [License text](http://www.boost.org/LICENSE_1_0.txt)
* [Discussion](http://www.boost.org/users/license.html)

# About
MathXK is maintained by Rocco Romeo. 

To report bugs or request new features create a new issue. Contributions are welcome.

For other inquiries, contact Rocco Romeo: rromeo@live.ca

# Acknowledgements
Thanks to the people at Boost Math who maintain the terrific [Boost Math](http://boost.org/libs/math) library. (Special thanks to John Maddock & Paul Bristow!)




